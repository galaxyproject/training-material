
> <solution-title>``client/src/entry/analysis/router.js``</solution-title>
> 
> Possible changes to file ``client/src/entry/analysis/router.js``:
> 
> ```diff
> index e4b3ce87cc..73332bf4ac 100644
> --- a/client/src/entry/analysis/router.js
> +++ b/client/src/entry/analysis/router.js
> @@ -22,6 +22,7 @@ import Grid from "components/Grid/Grid";
>  import GridShared from "components/Grid/GridShared";
>  import GridHistory from "components/Grid/GridHistory";
>  import HistoryImport from "components/HistoryImport";
> +import { FavoriteExtensions } from "components/User/FavoriteExtensions/index";
>  import HistoryView from "components/HistoryView";
>  import InteractiveTools from "components/InteractiveTools/InteractiveTools";
>  import InvocationReport from "components/Workflow/InvocationReport";
> @@ -299,6 +300,11 @@ export function getRouter(Galaxy) {
>                          props: true,
>                          redirect: redirectAnon(),
>                      },
> +                    {
> +                        path: "user/favorite/extensions",
> +                        component: FavoriteExtensions,
> +                        redirect: redirectAnon(),
> +                    },
>                      {
>                          path: "visualizations",
>                          component: VisualizationsList,
> ```
{: .solution }
