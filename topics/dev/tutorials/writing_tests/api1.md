
> <solution-title></solution-title>
> 
> Example of a basic API test.
> 
> ```python
> from ._framework import ApiTestCase
> 
> 
> class MyTutorialApiTestCase(ApiTestCase):
> 
>     def test_version_is_current(self):
>         response = self._get("version")
>         response.raise_for_status()
>         data = response.json()
>         assert data["version_major"] == "22.09"  # Assuming you are on the "dev" branch
> ```
{: .solution }
