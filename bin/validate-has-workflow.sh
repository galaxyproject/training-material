#!/bin/bash
# First find the disabled tutorials, we don't care if boxes are broken there
# Then find any 'raw' blockquote (missing any classes)
# and diff those, removing everything from the disabled list.
warning() {
	(>&2 echo "$(tput setab 2)$*$(tput sgr0)")
}

DISABLED_TUTOS="$(grep 'enable:.*false' -R topics/*/tutorials/*/tutorial.md -l | sed 's|/tutorial.md||')"


does_not_need_workflows() {
	cat | grep -v topics/dev | grep -v topics/admin | grep -v topics/contributing | grep -v  galaxy-interface | grep -v instructors
}


ALL_TUTOS="$(find topics/*/tutorials/ -maxdepth 1 -mindepth 1 | sort -u | does_not_need_workflows)"
HAS_WORKFLOWS="$(find topics/ | grep workflows | grep -v 'workflows$' | xargs -I{} dirname {} | xargs -I{} dirname {} | sort -u | does_not_need_workflows)"


echo "Missing Workflows"
diff <(echo "$ALL_TUTOS") <(echo "$HAS_WORKFLOWS") | grep '<' | grep -v "$DISABLED_TUTOS"
