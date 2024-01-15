#!/bin/bash

MATRIX_SERVER=${MATRIX_SERVER:-"https://matrix.org"}
ROOM_ID=${ROOM_ID:-'!TJRLNvfcbWbSRoUNpl:matrix.org'}  ## GTN Single Cell Maintainers
WANTED_TAGS=${WANTED_TAGS:-"scrna scrna-seq"}
MAX_REPLIES=${MAX_REPLIES:-1}
HTML_TYPE=${HTML_TYPE:-"bullets"}  ## "table"

## Result filters
OPTS=${OPTS:-"?ascending=true&order=activity"}

if [ -z "$MATRIX_ACCESS_TOKEN" ]; then
    echo "
This is a Matrix bot that scrapes Galaxy Help for certain tags and posts to
a Room for topics that have less than X replies. Run this maybe once a month.

Example Usage:

   MATRIX_ACCESS_TOKEN='123_123_123' \\
     MATRIX_SERVER='https://matrix.org' \\
     ROOM_ID='!123_132_123:matrix.org' \\
     WANTED_TAGS='tag1 tag2' \\
     MAX_REPLIES=1 \\
     HTML_TYPE='bullets' \\
     bash $0

Where:
  MATRIX_ACCESS_TOKEN  Can be found in your Matrix profile under
                       'All settings' -> 'Help & About' -> 'Access Token'

  MATRIX_SERVER        The name or base address of the Matrix server to
                       post to. Default is '$MATRIX_SERVER'

  ROOM_ID              The Room ID can be found in the URL of the room
                       usually following format '!123123123:matrix.org'.
                       Default is '$ROOM_ID'
                       NOTE: Single quotes are very important here.

  WANTED_TAGS          A space separated list of valid tags to find posts
                       at https://help.galaxyproject.org/
                       Default is \"$WANTED_TAGS\"

  MAX_REPLIES          Filter for posts that have less than or equal to
                       this many replies. Default is \"$MAX_REPLIES\"

  HTML_TYPE            Render either a 'table' or 'bullets'. HTML tables
                       look great in the browser but don't render well on
                       mobile. Default is \"$HTML_TYPE\"

  OPTS                 Extra arguments to append to help.galaxyproject.org
                       URL. Default is \"$OPTS\"
   " >&2
    exit 255
fi

function tag_to_tsv {
    ## For a given TAG, fetch from the help forum, extract and parse
    ## the table and produce a 4-column TSV output of Link, Title,
    ## Replies, Views.
    ##
    ## TODO: Add date too?
    local tag="$1"
    curl -s "https://help.galaxyproject.org/tag/${tag}${OPTS}" \
        | xmllint --noblanks --html --xpath '//tr[@class="topic-list-item"]/td/span' - 2>/dev/null \
        | sed -r 's|<span class[^>]+>||; s|</span>||; s|\s*<a.*href=\"([^\"]*)\" [^>]+>([^<]+)<.*|_ROW_\1\t\2|' \
        | tr '\n' '\t' \
        | sed -r 's|\s\s\s*|\t|g; s|_ROW_|\n|g'
}

function alltags_to_tsv {
    ## For all wanted tags, populate a 4-column TSV output of Link,
    ## Title, Replies, Views, and return the path of the table.
    local fetch_tags=$WANTED_TAGS
    local tmp_tsv
    tmp_tsv=$(mktemp --suffix=".tsv")
    for tag in ${fetch_tags}; do
        tag_to_tsv "$tag" >> "$tmp_tsv";
    done
    ## No duplicates, no blanks, no duplicate delimiters,
    ## and sort by ascending reply count
    grep -v "^\s*$" "${tmp_tsv}" | sed 's|\t\t|\t|g' \
        | sort | uniq | sort -t $'\t' -nk 3 > "${tmp_tsv}".temp
    echo -e "Link\tTitle\tReplies\tViews" > "${tmp_tsv}"
    cat "${tmp_tsv}".temp  >> "${tmp_tsv}"
    rm "${tmp_tsv}".temp
    echo "${tmp_tsv}"
}

function filter_tsv {
    ## Filter a TSV file for maximum replies and then return the path
    ## of the new filtered table
    local tsv="$1"
    local tmp_tsv
    tmp_tsv=$(mktemp --suffix=".tsv")
    awk -F$'\t' -v replies="$MAX_REPLIES" '$3 <= replies' "${tsv}" > "${tmp_tsv}"
    echo "${tmp_tsv}"
}

function tsv_to_html {
    ## Convert a TSV table into HTML text that can be fed into a JSON
    local tsv="$1"
    if [ "$HTML_TYPE" = "table" ]; then
        awk -F$'\t' -v subtitle="Recent posts matching: <b>${WANTED_TAGS}</b>, with replies &le; ${MAX_REPLIES}" '\
BEGIN { print "<h1>Updates from Galaxy Help</h1>"subtitle"\n<table>\n<thead><tr><th>Topic</th><th>Replies</th><th>Views</th></tr></thead>\n<tbody>"} \
END { print "</tbody>\n</table>"} \
NR > 0 {print "<tr><td><a href=\""$1"\">"$2"</a></td><td>"$3"</td><td>"$4"</td></tr>"}' \
            "${tsv}" | tr '\n' ' ' | sed 's|"|\\"|g'
    else  ## bullets
        awk -F$'\t' -v subtitle="Recent posts matching: <b>${WANTED_TAGS}</b>, with replies &le; ${MAX_REPLIES}" '\
BEGIN { print "<h1>Updates from Galaxy Help</h1><br/><p>"subtitle"</p><ol>\n"} \
END { print "\n</ol>"} \
NR > 0 {print "<li><a href=\""$1"\">"$2"</a><ul><li>Replies: "$3" and Views: "$4"</li></ul></li>"}' \
            "${tsv}" | tr '\n' ' ' | sed 's|"|\\"|g'
    fi
}

function tsv_to_markdown {
    ## Convert a TSV table into Markdown text that can be fed into a JSON
    local tsv="$1"
    awk -F$'\t' -v subtitle="Recent posts matching: **${WANTED_TAGS}**, with replies ≤ ${MAX_REPLIES}" '\
BEGIN { print "## Updates from Galaxy Help\\n***"subtitle"***\\n"} \
NR > 0 {print "* ["$2"]("$1")\\n   * "$3" replies and "$4" views\\n"}' \
        "${tsv}" | tr '\n' ' ' | sed 's|"|\\"|g'
}

function md_and_html_to_json {
    ## Stuff the Markdown and HTML text content into a JSON.
    local md_text="$1"
    local html_text="$2"
    local tmp_json;
    tmp_json=$(mktemp --suffix=".json")
    ## See: https://spec.matrix.org/legacy/r0.0.0/client_server.html
    ##      "put-matrix-client-r0-rooms-roomid-send-eventtype-txnid"
    echo "{\
\"msgtype\":\"m.notice\", \
\"format\":\"org.matrix.custom.html\", \
\"body\": \"${md_text}\", \
\"formatted_body\": \"${html_text}\"}" > "${tmp_json}"
    echo "${tmp_json}"
}

function post_json_to_matrix {
    local json_file="$1"
    local txnid post_url
    txnid=$(date "+%Y%m%d%H%M${RANDOM:1:3}")   ## date-specific transaction ID
    MATRIX_SERVER=${MATRIX_SERVER%/}           ## remove trailing slash, if any
    ## Build curl
    post_url="${MATRIX_SERVER}/_matrix/client/r0/rooms/"
    post_url="${post_url}"${ROOM_ID}"/send/m.room.message/${txnid}"
    post_url="${post_url}?access_token=${MATRIX_ACCESS_TOKEN}"
    ## DEBUG:
    ## - curl "$post_url" -X PUT --data '{"msgtype":"m.text","body":"hello"}'
    curl "$post_url" -X PUT --data "$(cat ${json_file})"
}

function sanity_check {
    ## Assert that required binaries are in PATH
    local required_progs=( cat curl xmllint awk sed grep tr jq )
    local miss=""
    for prog in "${required_progs[@]}"; do
        if ! which "${prog}" 2>/dev/null >&2; then
            miss="$miss $prog"
        fi
    done
    if [ "$miss" != "" ]; then
        echo "Cannot run without:$miss"
        exit 255
    fi
}

## MAIN ##
sanity_check

main_tsv=$(filter_tsv "$(alltags_to_tsv)" )
if [[ $(wc -l < "${main_tsv}") == 0 ]]; then
    echo "Nothing new to post, aborting."  >&2
    exit 0
fi

main_mdwn_text=$(tsv_to_markdown "${main_tsv}")
main_html_text=$(tsv_to_html "${main_tsv}")

main_json_file=$(md_and_html_to_json "${main_mdwn_text}" "${main_html_text}")
if ! jq < "${main_json_file}" 2> /dev/null >&2; then
    echo "This is not a valid JSON, aborting." >&2
    echo "See: ${main_json_file}"  >&2
    exit 255
fi

post_json_to_matrix "${main_json_file}"
