
> ### {% icon solution %} ``client/src/entry/analysis/AnalysisRouter.js``
> 
> Possible changes to file ``client/src/entry/analysis/AnalysisRouter.js``:
> 
> ```diff
> index b409c8374f..61b01f85fa 100644
> --- a/client/src/entry/analysis/AnalysisRouter.js
> +++ b/client/src/entry/analysis/AnalysisRouter.js
> @@ -32,6 +32,7 @@ import InteractiveTools from "components/InteractiveTools/InteractiveTools.vue";
>  import WorkflowList from "components/Workflow/WorkflowList.vue";
>  import HistoryImport from "components/HistoryImport.vue";
>  import { HistoryExport } from "components/HistoryExport/index";
> +import { FavoriteExtensions } from "components/User/FavoriteExtensions/index";
>  import HistoryView from "components/HistoryView.vue";
>  import WorkflowInvocationReport from "components/Workflow/InvocationReport.vue";
>  import WorkflowRun from "components/Workflow/Run/WorkflowRun.vue";
> @@ -67,6 +68,7 @@ export const getAnalysisRouter = (Galaxy) => {
>              "(/)user(/)cloud_auth": "show_cloud_auth",
>              "(/)user(/)external_ids": "show_external_ids",
>              "(/)user(/)(:form_id)": "show_user_form",
> +            "(/)user(/)favorite(/)extensions": "show_user_favorite_extensions",
>              "(/)pages(/)create(/)": "show_pages_create",
>              "(/)pages(/)edit(/)": "show_pages_edit",
>              "(/)pages(/)sharing(/)": "show_pages_sharing",
> @@ -157,6 +159,10 @@ export const getAnalysisRouter = (Galaxy) => {
>              this.page.display(new FormWrapper.View(_.extend(model[form_id], { active_tab: "user" })));
>          },
>  
> +        show_user_favorite_extensions: function () {
> +            this._display_vue_helper(FavoriteExtensions, {});
> +        },
> +
>          show_interactivetool_list: function () {
>              this._display_vue_helper(InteractiveTools);
>          },
> ```
{: .solution }
