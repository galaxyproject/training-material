
> <solution-title>``client/src/components/User/FavoriteExtensions/List.test.js``</solution-title>
> 
> Possible changes to file ``client/src/components/User/FavoriteExtensions/List.test.js``:
> 
> ```diff
> new file mode 100644
> index 0000000000..7efa1c76ec
> --- /dev/null
> +++ b/client/src/components/User/FavoriteExtensions/List.test.js
> @@ -0,0 +1,97 @@
> +import axios from "axios";
> +import MockAdapter from "axios-mock-adapter";
> +import { mount } from "@vue/test-utils";
> +import { getLocalVue } from "jest/helpers";
> +import List from "./List";
> +import flushPromises from "flush-promises";
> +
> +jest.mock("app");
> +
> +const localVue = getLocalVue();
> +const propsData = {};
> +
> +describe("User/FavoriteExtensions/List.vue", () => {
> +    let wrapper;
> +    let axiosMock;
> +
> +    beforeEach(() => {
> +        axiosMock = new MockAdapter(axios);
> +        axiosMock.onGet("api/datatypes").reply(200, ["fasta", "fastq"]);
> +        axiosMock.onGet("api/users/current/favorites/extensions").reply(200, ["fasta"]);
> +    });
> +
> +    afterEach(() => {
> +        axiosMock.restore();
> +    });
> +
> +    it("should start loading and then render list", async () => {
> +        wrapper = mount(List, {
> +            propsData,
> +            localVue,
> +        });
> +        expect(wrapper.vm.loading).toBeTruthy();
> +        await flushPromises();
> +        expect(wrapper.vm.loading).toBeFalsy();
> +        const el = wrapper.find("#favorite-extensions");
> +        expect(el.exists()).toBe(true);
> +        expect(el.find("[data-extension-target]").exists()).toBe(true);
> +    });
> +
> +    it("should mark favorite and not favorite with different links", async () => {
> +        wrapper = mount(List, {
> +            propsData,
> +            localVue,
> +        });
> +        await flushPromises();
> +        const el = wrapper.find("#favorite-extensions");
> +        const els = el.findAll("[data-extension-target]");
> +        expect(els.length).toBe(2);
> +        const fastaEntry = els.at(0);
> +        expect(fastaEntry.attributes("data-extension-target")).toBe("fasta");
> +        expect(fastaEntry.find(".unmark-favorite").exists()).toBe(true);
> +        expect(fastaEntry.find(".mark-favorite").exists()).toBe(false);
> +
> +        const fastqEntry = els.at(1);
> +        expect(fastqEntry.attributes("data-extension-target")).toBe("fastq");
> +        expect(fastqEntry.find(".mark-favorite").exists()).toBe(true);
> +        expect(fastqEntry.find(".unmark-favorite").exists()).toBe(false);
> +    });
> +
> +    it("should post to mark favorites", async () => {
> +        wrapper = mount(List, {
> +            propsData,
> +            localVue,
> +        });
> +        await flushPromises();
> +        const el = wrapper.find("#favorite-extensions");
> +        const els = el.findAll("[data-extension-target]");
> +        const fastqEntry = els.at(1);
> +        const markFastq = fastqEntry.find(".mark-favorite a");
> +        expect(markFastq.exists()).toBe(true);
> +
> +        axiosMock.onPost("api/users/current/favorites/extensions/fastq").reply(200, "fastq");
> +        axiosMock.onGet("api/users/current/favorites/extensions").reply(200, ["fasta", "fastq"]);
> +        markFastq.trigger("click");
> +        await flushPromises();
> +        expect(wrapper.vm.favoriteExtensions.indexOf("fastq") >= 0).toBe(true);
> +    });
> +
> +    it("should delete to unmark favorites", async () => {
> +        wrapper = mount(List, {
> +            propsData,
> +            localVue,
> +        });
> +        await flushPromises();
> +        const el = wrapper.find("#favorite-extensions");
> +        const els = el.findAll("[data-extension-target]");
> +        const fastaEntry = els.at(0);
> +        const unmarkFasta = fastaEntry.find(".unmark-favorite a");
> +        expect(unmarkFasta.exists()).toBe(true);
> +
> +        axiosMock.onDelete("api/users/current/favorites/extensions/fasta").reply(200);
> +        axiosMock.onGet("api/users/current/favorites/extensions").reply(200, []);
> +        unmarkFasta.trigger("click");
> +        await flushPromises();
> +        expect(wrapper.vm.favoriteExtensions.indexOf("fasta") < 0).toBe(true);
> +    });
> +});
> ```
{: .solution }
