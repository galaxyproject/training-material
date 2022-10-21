
> <solution-title>``client/src/components/User/FavoriteExtensions/List.vue``</solution-title>
> 
> Possible changes to file ``client/src/components/User/FavoriteExtensions/List.vue``:
> 
> ```diff
> new file mode 100644
> index 0000000000..053b354f63
> --- /dev/null
> +++ b/client/src/components/User/FavoriteExtensions/List.vue
> @@ -0,0 +1,78 @@
> +<template>
> +    <div class="favorite-extensions-card">
> +        <b-alert variant="error" show v-if="errorMessage"> </b-alert>
> +        <loading-span v-if="loading" message="Loading favorite extensions" />
> +        <ul id="favorite-extensions" v-else>
> +            <li v-for="extension in extensions" :key="extension" :data-extension-target="extension">
> +                <span
> +                    class="favorite-link unmark-favorite"
> +                    v-if="favoriteExtensions.indexOf(extension) >= 0"
> +                    title="Unmark as favorite">
> +                    <a href="#" @click="unmarkAsFavorite(extension)">(X)</a>
> +                </span>
> +                <span class="favorite-link mark-favorite" v-else title="Mark as favorite">
> +                    <a href="#" @click="markAsFavorite(extension)">(+)</a>
> +                </span>
> +            </li>
> +        </ul>
> +    </div>
> +</template>
> +
> +<script>
> +import axios from "axios";
> +import { getGalaxyInstance } from "app";
> +import LoadingSpan from "components/LoadingSpan";
> +import { errorMessageAsString } from "utils/simple-error";
> +
> +export default {
> +    components: {
> +        LoadingSpan,
> +    },
> +    data() {
> +        const Galaxy = getGalaxyInstance();
> +        return {
> +            datatypesUrl: `${Galaxy.root}api/datatypes`,
> +            favoriteExtensionsUrl: `${Galaxy.root}api/users/current/favorites/extensions`,
> +            extensions: null,
> +            favoriteExtensions: null,
> +            errorMessage: null,
> +        };
> +    },
> +    created() {
> +        this.loadDatatypes();
> +        this.loadFavorites();
> +    },
> +    computed: {
> +        loading() {
> +            return this.extensions == null || this.favoriteExtensions == null;
> +        },
> +    },
> +    methods: {
> +        loadDatatypes() {
> +            axios
> +                .get(this.datatypesUrl)
> +                .then((response) => {
> +                    this.extensions = response.data;
> +                })
> +                .catch(this.handleError);
> +        },
> +        loadFavorites() {
> +            axios
> +                .get(this.favoriteExtensionsUrl)
> +                .then((response) => {
> +                    this.favoriteExtensions = response.data;
> +                })
> +                .catch(this.handleError);
> +        },
> +        markAsFavorite(extension) {
> +            axios.post(`${this.favoriteExtensionsUrl}/${extension}`).then(this.loadFavorites).catch(this.handleError);
> +        },
> +        unmarkAsFavorite(extension) {
> +            axios.delete(`${this.favoriteExtensionsUrl}/${extension}`).then(this.loadFavorites).catch(this.handleError);
> +        },
> +        handleError(error) {
> +            this.errorMessage = errorMessageAsString(error);
> +        },
> +    },
> +};
> +</script>
> ```
{: .solution }
