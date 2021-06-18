
> ### {% icon solution %} ``client/src/layout/menu.js``
> 
> Possible changes to file ``client/src/layout/menu.js``:
> 
> ```diff
> index 062191b492..de91d9670c 100644
> --- a/client/src/layout/menu.js
> +++ b/client/src/layout/menu.js
> @@ -314,6 +314,11 @@ export function fetchMenu(options = {}) {
>                      url: "workflows/invocations",
>                      target: "__use_router__",
>                  },
> +                {
> +                    title: _l("Favorite Extensions"),
> +                    url: "/user/favorite/extensions",
> +                    target: "__use_router__",
> +                },
>              ],
>          };
>          if (Galaxy.config.visualizations_visible) {
> ```
{: .solution }
