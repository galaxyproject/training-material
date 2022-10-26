
> <solution-title>``client/src/layout/menu.js``</solution-title>
> 
> Possible changes to file ``client/src/layout/menu.js``:
> 
> ```diff
> index b4f6c46a2f..d2eab16f6f 100644
> --- a/client/src/layout/menu.js
> +++ b/client/src/layout/menu.js
> @@ -303,6 +303,11 @@ export function fetchMenu(options = {}) {
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
