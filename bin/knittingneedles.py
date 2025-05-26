import copy
import re

BOXES = r"^([\s>]*>[\s>]*)"
BOX_PREFIX = r"\s*{% raw %}"
BOX_OPEN = r"\s*```diff"
BOX_CLOSE = r'\s*{: data-commit="([^"]*)"([^}]*)}'
BOX_CLOSE_ALL = r'\s*{:\s*(.spoken|data-commit|.code-in\s*data-cmd).*}'
WHITESPACE = r"^(\s*)"

CMD_OPEN = r"\s*```bash"
CMD_CLOSE = r'\s*{: data-cmd="true"([^}]*)}'

TEST_OPEN = r"\s*```bash"
TEST_CLOSE = r'\s*{: data-test="true"([^}]*)}'


def stripN(line, count):
    c = copy.copy(line)
    for i in range(count):
        c = c.lstrip()
        c = c.lstrip(">")
    removed = len(line) - len(c)
    return (c, line[0:removed])


def removeWhitespacePrefix(lines):
    # Obtain whitespace amount
    whitespace_prefix = [len(re.match(WHITESPACE, line).group(1)) for line in lines]
    # Remove zeros (blank lines have no whitespace, not even the standard of
    # the rest)
    whitespace_prefix = [x for x in whitespace_prefix if x != 0]
    amount = min(whitespace_prefix)
    diff = [x[amount:] for x in lines]
    return (amount, diff)


def extractCommitMsg(lines):
    return re.match(BOX_CLOSE, lines[-1].strip()).group(1)
