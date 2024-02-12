curl https://training.galaxyproject.org/sitemap.xml | xpath -q -e '//url/loc/text()' | \
	grep -v '/training-material/api/' | \
	egrep -v '\.ya?ml$' | \
	grep -v 'slides-plain.html$' | \
	sed 's|https://training.galaxyproject.org|http://localhost:4002|' > /tmp/urls
for x in `cat /tmp/urls`; do 
	resp=$(curl -s -w "%{http_code}\n" $x -o /dev/null);
	if (( resp != 200 )); then
		echo $x; 
	fi;
done
