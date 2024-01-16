
> <solution-title></solution-title>
> 
> Alternative test for creating an object. Note the `setUp` method and the new import.
> 
> ```python
> from galaxy_test.base.populators import DatasetPopulator
> from ._framework import ApiTestCase
> 
> 
> class MyTutorialApiTestCase(ApiTestCase):
>
>     def setUp(self):
>         super().setUp()
>         self.dataset_populator = DatasetPopulator(self.galaxy_interactor)
> 
> # [code omitted]
> 
>     def test_create_role2(self):
>         # prepare new role
>         name = self.dataset_populator.get_random_name()
>         description = 'description of this cool role'
>         payload = {
>             "name": name,
>             "description": description,
>         }
> 
>         # verify new role does not exist
>         response = self._get("roles")
>         response.raise_for_status()
>         data = response.json()
>         assert not any(role["name"] == name for role in data)
> 
>         # add new role
>         response = self._post("roles", payload, admin=True, json=True)
>         response.raise_for_status()
>         role = response.json()
>         assert role["name"] == name
>         assert role["description"] == description
> 
>         # verify role has been added
>         response = self._get("roles")
>         response.raise_for_status()
>         data = response.json()
>         assert any(role["name"] == name for role in data)
> ```
{: .solution }
