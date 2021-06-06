
> ### {% icon solution %} ``lib/galaxy/webapps/galaxy/buildapp.py``
> 
> Possible changes to file ``lib/galaxy/webapps/galaxy/buildapp.py``:
> 
> ```diff
> index 7fbe248c66..76fe7cd819 100644
> --- a/lib/galaxy/webapps/galaxy/buildapp.py
> +++ b/lib/galaxy/webapps/galaxy/buildapp.py
> @@ -155,6 +155,7 @@ def app_pair(global_conf, load_app_kwds=None, wsgi_preflight=True, **kwargs):
>      webapp.add_client_route('/tours/{tour_id}')
>      webapp.add_client_route('/user')
>      webapp.add_client_route('/user/{form_id}')
> +    webapp.add_client_route('/user/favorite/extensions')
>      webapp.add_client_route('/visualizations')
>      webapp.add_client_route('/visualizations/edit')
>      webapp.add_client_route('/visualizations/sharing')
> ```
{: .solution }
