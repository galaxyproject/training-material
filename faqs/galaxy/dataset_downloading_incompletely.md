| title   |  area  | layout |	box_type | contributors|
|---------|:------|:-------|:------------|:------------|
|Dataset downloading incompletely?| Datasets |faq| tip| Andreas_sune|



1. Try using the [Google Chrome web browser](https://www.google.com/chrome/?brand=BNSD&gclid=Cj0KCQjwxtSSBhDYARIsAEn0thSKnKsxxKVx6hXa7XtuUD7WrhkA0BWbQE3Q2a47Y5V2S69Sh2qgqowaAtqfEALw_wcB&gclsrc=aw.ds). Sometimes Chrome itself better supports continuous data transfers. 


2.  Use the command-line option to download your dataset:

	* Copy the link location of your dataset.
	* Open the terminal window of your computer.
	*  Use either  [wget](https://www.gnu.org/software/wget/manual/html_node/Download-Options.html#Download-Options) command :
		```
		$ wget -O '<link>'
		$ wget -O --no-check-certificate '<link>' # ignore SSL certificate warnings
		$ wget -c '<link>'  # continue an interrupted download
	   ```
	* Or [curl](https://en.wikipedia.org/wiki/CURL) command :
	    ```
	    $ curl -o outfile '<link>' 
	    $ curl -o outfile --insecure '<link>'     # ignore SSL certificate warnings
        $ curl -C - -o outfile '<link>'           # continue an interrupted download
        ```
	

**Note:**
* The commands have many options. These are examples commonly used with Galaxy.
* "$" indicates the terminal prompt.
			