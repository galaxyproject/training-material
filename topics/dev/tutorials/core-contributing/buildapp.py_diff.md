
> <solution-title>``lib/galaxy/webapps/galaxy/buildapp.py``</solution-title>
> 
> Possible changes to file ``lib/galaxy/webapps/galaxy/buildapp.py``:
> 
> ```diff
> index 83757e5307..c9e0feeb77 100644
> --- a/lib/galaxy/webapps/galaxy/buildapp.py
> +++ b/lib/galaxy/webapps/galaxy/buildapp.py
> @@ -213,6 +213,7 @@ def app_pair(global_conf, load_app_kwds=None, wsgi_preflight=True, **kwargs):
>      webapp.add_client_route("/tours/{tour_id}")
>      webapp.add_client_route("/user")
>      webapp.add_client_route("/user/{form_id}")
> +    webapp.add_client_route('/user/favorite/extensions')
>      webapp.add_client_route("/welcome/new")
>      webapp.add_client_route("/visualizations")
>      webapp.add_client_route("/visualizations/edit")
> -- 
> 2.30.1 (Apple Git-130)
> 
> ```
{: .solution }
