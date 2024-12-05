#!/usr/bin/env python3
"""Simple curses-based text editor.

This code is based on https://github.com/tdryer/editor by tdryer. The original code is licensed under MIT:

The MIT License (MIT)

Copyright (c) 2014 Tom Dryer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import time
import argparse
from contextlib import contextmanager
import whatthepatch
import sys
import curses
import os
CHAR_ESC = 27
CHAR_BKSP = 127


class Buffer(object):
    """The basic data structure for editable text.

    The buffer is column and row oriented. Column and row numbers start with 0.
    A buffer always has at least one row. All positions within a buffer specify
    a position between characters.
    """

    def __init__(self, text=''):
        """Create a new Buffer, optionally initialized with text."""
        self._lines = text.split('\n')

    def get_lines(self):
        """Return list of lines in the buffer."""
        return list(self._lines) # return a copy

    def _check_point(self, row, col):
        """Raise ValueError if the given row and col are not a valid point."""
        if row < 0 or row > len(self._lines) - 1:
            raise ValueError("Invalid row: '{}'".format(row))
        cur_row = self._lines[row]
        if col < 0 or col > len(cur_row):
            raise ValueError("Invalid col: '{}'".format(col))

    def set_text(self, row1, col1, row2, col2, text):
        """Set the text in the given range.

        The end of the range is exclusive (to allow inserting text without
        removing a single character). Column numbers are positions between
        characters.

        Raises ValueError if the range is invalid.
        """
        # TODO check that point2 is after or the same as point1
        self._check_point(row1, col1)
        self._check_point(row2, col2)

        line = self._lines[row1][:col1] + text + self._lines[row2][col2:]
        self._lines[row1:row2+1] = line.split('\n')


class EditorGUI(object):

    def __init__(self, stdscr, filename, stream=None, speed=0, nosave=False, seshlen=0):
        """Create the GUI with curses screen and optional filename to load."""
        self._stdscr = stdscr

        # if filename already exists, try to load from it
        text = ''
        if filename != None and os.path.isfile(filename):
            with open(filename) as f:
                text = f.read()

        self._filename = filename
        self._buf = Buffer(text)
        self._row = 0
        self._col = 0
        self._scroll_top = 0 # the first line number in the window
        self._mode = 'normal'
        self._message = ''
        self._will_exit = False
        # Extras
        self._char_stream = stream
        self._speed = speed
        self._nosave = nosave
        self._seshlen = seshlen
        self._start_time = time.time()

    def _draw_gutter(self, num_start, num_rows, last_line_num):
        """Draw the gutter, and return the gutter width."""
        line_nums = list(range(num_start, num_start + num_rows))
        assert len(line_nums) == num_rows
        gutter_width = max(3, len(str(last_line_num))) + 1
        for y, line_num in enumerate(line_nums):
            if line_num > last_line_num:
                text = '~'.ljust(gutter_width)
            else:
                text = '{} '.format(line_num).rjust(gutter_width)
            self._stdscr.addstr(y, 0, text, curses.A_REVERSE)
        return gutter_width

    def _draw(self):
        """Draw the GUI."""
        self._stdscr.erase()
        height = self._stdscr.getmaxyx()[0]
        width = self._stdscr.getmaxyx()[1]
        # self._draw_status_line(0, height - 1, width)
        self._draw_text(0, 0, width, height - 1)
        self._stdscr.refresh()

    def _draw_status_line(self, left, top, width):
        """Draw the status line."""
        # TODO: can't write to bottom right cell
        mode = '{} {} {}'.format(self._filename, self._mode.upper(),
                                 self._message).ljust(width - 1)
        self._stdscr.addstr(top, left, mode, curses.A_REVERSE)
        position = 'LN {}:{} '.format(self._row + 1, self._col + 1)
        self._stdscr.addstr(top, left + width - 1 - len(position), position,
                            curses.A_REVERSE)

    def _get_num_wrapped_lines(self, line_num, width):
        """Return the number of lines the given line number wraps to."""
        return len(self._get_wrapped_lines(line_num, width))

    def _get_wrapped_lines(self, line_num, width, convert_nonprinting=True):
        """Return the wrapped lines for the given line number."""
        def wrap_text(text, width):
            """Wrap string text into list of strings."""
            if text == '':
                yield ''
            else:
                for i in range(0, len(text), width):
                    yield text[i:i + width]
        assert line_num >= 0, 'line_num must be > 0'
        line = self._buf.get_lines()[line_num]
        if convert_nonprinting:
            line = self._convert_nonprinting(line)
        return list(wrap_text(line, width))

    def _scroll_bottom_to_top(self, bottom, width, height):
        """Return the first visible line's number so bottom line is visible."""
        def verify(top):
            """Verify the result of the parent function is correct."""
            rows = [list(self._get_wrapped_lines(n, width))
                    for n in range(top, bottom + 1)]
            num_rows = sum(len(r) for r in rows)
            assert top <= bottom, ('top line {} may not be below bottom {}'
                                   .format(top, bottom))
            assert num_rows <= height, (
                '{} rows between {} and {}, but only {} remaining. rows are {}'
                .format(num_rows, top, bottom, height, rows))

        top, next_top = bottom, bottom
        # distance in number of lines between top and bottom
        distance = self._get_num_wrapped_lines(bottom, width)

        # move top upwards as far as possible
        while next_top >= 0 and distance <= height:
            top = next_top
            next_top -= 1
            distance += self._get_num_wrapped_lines(max(0, next_top), width)

        verify(top)
        return top

    def _scroll_to(self, line_num, width, row_height):
        """Scroll so the line with the given number is visible."""
        # lowest scroll top that would still keep line_num visible
        lowest_top = self._scroll_bottom_to_top(line_num, width, row_height)

        if line_num < self._scroll_top:
            # scroll up until line_num is visible
            self._scroll_top = line_num
        elif self._scroll_top < lowest_top:
            # scroll down to until line_num is visible
            self._scroll_top = lowest_top

    @staticmethod
    def _convert_nonprinting(text):
        """Replace nonprinting character in text."""
        # TODO: it would be nice if these could be highlighted when displayed
        res = []
        for char in text:
            i = ord(char)
            if char == '\t':
                res.append('->  ')
            elif i < 32 or i > 126:
                res.append('<{}>'.format(hex(i)[2:]))
            else:
                res.append(char)
        return ''.join(res)

    def _draw_text(self, left, top, width, height):
        """Draw the text area."""
        # TODO: handle single lines that occupy the entire window
        highest_line_num = len(self._buf.get_lines())
        gutter_width = max(3, len(str(highest_line_num))) + 1
        line_width = width - gutter_width # width to which text is wrapped
        cursor_y, cursor_x = None, None # where the cursor will be drawn

        # set scroll_top so the cursor is visible
        self._scroll_to(self._row, line_width, height)

        line_nums = list(range(self._scroll_top, highest_line_num))
        cur_y = top
        trailing_char = '~'

        for line_num in line_nums:

            # if there are no more rows left, break
            num_remaining_rows = top + height - cur_y
            if num_remaining_rows == 0:
                break

            # if all the wrapped lines can't fit on screen, break
            wrapped_lines = self._get_wrapped_lines(line_num, line_width)
            if len(wrapped_lines) > num_remaining_rows:
                trailing_char = '@'
                break

            # calculate cursor position if cursor must be on this line
            if line_num == self._row:
                lines = self._get_wrapped_lines(line_num, line_width,
                                                convert_nonprinting=False)
                real_col = len(self._convert_nonprinting(
                    ''.join(lines)[:self._col])
                )
                cursor_y = cur_y + real_col / line_width
                cursor_x = left + gutter_width + real_col % line_width

            # draw all the wrapped lines
            for n, wrapped_line in enumerate(wrapped_lines):
                if n == 0:
                    gutter = '{} '.format(line_num + 1).rjust(gutter_width)
                else:
                    gutter = ' ' * gutter_width
                self._stdscr.addstr(cur_y, left, gutter, curses.A_REVERSE)
                self._stdscr.addstr(cur_y, left + len(gutter), wrapped_line)
                cur_y += 1

        # draw empty lines
        for cur_y in range(cur_y, top + height):
            gutter = trailing_char.ljust(gutter_width)
            self._stdscr.addstr(cur_y, left, gutter)

        # position the cursor
        assert cursor_x != None and cursor_y != None
        self._stdscr.move(int(cursor_y) + 0, int(cursor_x) + 0)

    def _handle_normal_keypress(self, char):
        """Handle a keypress in normal mode."""
        if char == ord('q'): # quit
            self._will_exit = True
        elif char == ord('j'): # down
            self._row += 1
        elif char == ord('k'): # up
            self._row -= 1
        elif char == ord('h'): # left
            self._col -= 1
        elif char == ord('l'): # right
            self._col += 1
        elif char == ord('0'): # move to beginning of line
            self._col = 0
        elif char == ord('$'): # move to end of line
            cur_line_len = len(self._buf.get_lines()[self._row])
            self._col = cur_line_len - 1
        elif char == ord('x'): # delete a character
            self._buf.set_text(self._row, self._col, self._row,
                                self._col + 1, '')
        elif char == ord('i'): # enter insert mode
            self._mode = "insert"
        elif char == ord('a'): # enter insert mode after cursor
            self._mode = "insert"
            self._col += 1
        elif char == ord('o'): # insert line after current
            cur_line_len = len(self._buf.get_lines()[self._row])
            self._buf.set_text(self._row, cur_line_len, self._row,
                               cur_line_len, '\n')
            self._row += 1
            self._col = 0
            self._mode = "insert"
        elif char == ord('O'): # insert line before current
            self._buf.set_text(self._row, 0, self._row, 0, '\n')
            self._col = 0
            self._mode = "insert"
        elif char == ord('w'): # write file
            if self._filename == None:
                self._message = 'Can\'t write file without filename.'
            else:
                try:
                    with open(self._filename, 'w') as f:
                        f.write('\n'.join(self._buf.get_lines()))
                except IOError as e:
                    self._message = ('Failed to write file \'{}\': {}'
                                     .format(self._filename, e))
        else:
            self._message = 'Unknown key: {}'.format(char)

    def _handle_insert_keypress(self, char):
        """Handle a keypress in insert mode."""
        if char == CHAR_ESC:
            # leaving insert mode moves cursor left
            if self._mode == 'insert':
                self._col -= 1
            self._mode = "normal"
        elif char == CHAR_BKSP: # backspace
            if self._col == 0 and self._row == 0:
                pass # no effect
            elif self._col == 0:
                # join the current line with the previous one
                prev_line = self._buf.get_lines()[self._row - 1]
                cur_line = self._buf.get_lines()[self._row]
                self._buf.set_text(self._row - 1, 0, self._row,
                                    len(cur_line), prev_line + cur_line)
                self._col = len(prev_line)
                self._row -= 1
            else:
                # remove the previous character
                self._buf.set_text(self._row, self._col - 1, self._row,
                                    self._col, '')
                self._col -= 1
        else:
            self._message = ('inserted {} at row {} col {}'
                             .format(char, self._row, self._col))
            self._buf.set_text(self._row, self._col, self._row,
                                self._col, chr(char))
            if chr(char) == '\n':
                self._row += 1
                self._col = 0
            else:
                self._col += 1

    def main(self):
        """GUI main loop."""

        # Reverse to permit popping
        if self._char_stream is not None:
            stream = self._char_stream[::-1]
            automated = True
        else:
            automated = False

        while not self._will_exit:
            self._draw()
            self._message = ''

            time.sleep(self._speed)
            if automated:
                if len(stream) > 0:
                    if len(stream) == len(self._char_stream):
                        # Pause before starting
                        time.sleep(2)
                    char = stream.pop()
                else:
                    # Push write and quit
                    stream.append(ord('q'))
                    if not self._nosave:
                        stream.append(ord('w'))
                    # Pause before writing/exiting
                    now = time.time()
                    elapsed_time = now - self._start_time
                    if self._seshlen > elapsed_time:
                        time.sleep(self._seshlen - elapsed_time)
                    else:
                        time.sleep(2)
            else:
                char = self._stdscr.getch()

            if self._mode == 'normal':
                self._handle_normal_keypress(char)
            elif self._mode == 'insert':
                self._handle_insert_keypress(char)

            # TODO: get rid of this position clipping
            num_lines = len(self._buf.get_lines())
            # This is a workaround. If we press 'j' on the last line, it will
            # add another new line for us. This works around having to figure
            # out if we're on the last line to add a new line.
            row_max = min(num_lines - 1, max(0, self._row))
            if self._row > row_max:
                self._buf._lines.append('')
            self._row = min(num_lines, max(0, self._row))

            # on empty lines, still allow col 1
            num_cols = max(1, len(self._buf.get_lines()[self._row]))
            # in insert mode, allow using append after the last char
            if self._mode == 'insert':
                num_cols += 1
            self._col = min(num_cols - 1, max(0, self._col))


@contextmanager
def use_curses():
    """Context manager to set up and tear down curses."""
    stdscr = curses.initscr()
    curses.noecho() # do not echo keys
    curses.cbreak() # don't wait for enter
    try:
        yield stdscr
    finally:
        # clean up and exit
        curses.nocbreak()
        stdscr.keypad(0)
        curses.echo()
        curses.endwin()


def move(target, current):
    n = target - current
    if target > current:
        return (n, ['j'] * n)
    elif target < current:
        return (n, ['k'] * n)
    else:
        return (n, []) # On correct line


def curses_main():
    """Start the curses GUI."""
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--diff', type=argparse.FileType('r'), help="Path to the diff to apply")
    parser.add_argument('--file', type=argparse.FileType('r'), help="Path to file to edit (conflicts with --dif)")
    parser.add_argument('--speed', type=float, default=0.01, help="Time to sleep between button presses")
    parser.add_argument('--session-min-length', type=float, default=0, help="The minimum time of the recording session (for syncing with audio.) Will sleep until this time has been reached.")
    parser.add_argument('--debug', action='store_true', help="Print out character stream and exit")
    parser.add_argument('--nosave', action='store_true', help="Do not save the output")
    args = parser.parse_args()

    stream = None
    fn = None
    if args.file:
        fn = args.file.name

    if args.diff:
        p = list(whatthepatch.parse_patch(args.diff.read()))

        if len(p) != 1:
            raise Exception("Uhh can't parse this")

        (header, changes, text) = p[0]
        stream = []
        fn = header.new_path
        if header.new_path == '/dev/null':
            raise Exception('Removing files is unsupported!')
        elif header.old_path == '/dev/null' or header.old_path == header.new_path:
            # If there is a set of folders, we need to make it.
            if '/' in header.new_path:
                directory = os.path.dirname(header.new_path)
                os.makedirs(directory, exist_ok=True)

            current_line = 1
            line_delta = 0
            for c in changes:
                if c.new is None:
                    if args.debug:
                        stream.append(f'CASE A/- move {c.old + line_delta}<-{current_line}')
                    # print(f'removing line {c.old}: {c.line}')
                    (motion_count, motions) = move(c.old + line_delta, current_line)
                    stream.extend(motions)
                    current_line += motion_count

                    stream.extend(['x'] * len(c.line))
                    stream.extend(['i', chr(CHAR_BKSP), chr(CHAR_ESC), '0'])
                    current_line -= 1
                    line_delta -= 1

                elif c.old is None:
                    if args.debug:
                        stream.append(f'CASE B/+ move {c.new}<-{current_line}')
                    (motion_count, motions) = move(c.new, current_line)
                    stream.extend(motions)
                    current_line = c.new

                    stream.append('O') # Enter edit mode
                    stream.extend(c.line)
                    stream.append(chr(CHAR_ESC)) # Return to normal
                    stream.append('0') # Ensure at start of line
                else:
                    line_delta = c.new - c.old
                    # if args.debug:
                        # stream.append('CASE C/=')
                    # stream.extend(move(c.new, current_line))
                    # current_line = c.new

                    pass
                    # print(f'No change?? {c.old} {c.new} {c.line}')
                if args.debug:
                    print(c)
                    stream.append(f'DEBUG: {current_line}/{line_delta}')
                    print(stream)
                    stream = []
        else:
            raise Exception("Cannot handle renames")

        if args.debug:
            print(stream)
            sys.exit()

        # Just in case.
        stream.append(chr(CHAR_ESC))
        stream = list(map(ord, stream))

    with use_curses() as stdscr:
        gui = EditorGUI(stdscr, fn, stream=stream, speed=args.speed, nosave=args.nosave, seshlen=args.session_min_length)
        gui.main()


if __name__ == '__main__':
    curses_main()
