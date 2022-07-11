
> ### {% icon solution %} Solution
> 
> API test for verifying the list of datatypes.
> 
> ```python
> from ._framework import ApiTestCase
> 
> 
> NUMBER_OF_DATATYPES = 446
> 
> 
> class MyTutorialApiTestCase(ApiTestCase):
> 
>     def test_version_is_current(self):
>         response = self._get("version")
>         response.raise_for_status()
>         data = response.json()
>         assert data["version_major"] == "22.09"
> 
>     def test_list_datatypes(self):
>         response = self._get("datatypes")
>         response.raise_for_status()
>         data = response.json()
>         assert len(data) == NUMBER_OF_DATATYPES
>         assert 'fasta' in data
>         assert 'totally invalid' not in data
> ```
{: .solution }


