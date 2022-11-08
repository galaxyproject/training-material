
> <solution-title></solution-title>
> 
> Final code.
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
>     def test_version_is_current(self):
>         response = self._get("version")
>         response.raise_for_status()
>         data = response.json()
>         assert data["version_major"] == "22.09"
> 
>     def test_create_role(self):
>         # prepare new role
>         name = self.dataset_populator.get_random_name()
>         description = 'description of this cool role'
>         payload = {
>             "name": name,
>             "description": description,
>         }
>         # add new role and verify
>         response = self._post("roles", payload, admin=True, json=True)
>         response.raise_for_status()
>         role = response.json()
>         assert role["name"] == name
>         assert role["description"] == description
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
> 
>     def test_search_dataset_by_history(self):
>         self.history1_id = self.dataset_populator.new_history()
>         self.history2_id = self.dataset_populator.new_history()
> 
>         history1_datasets, history2_datasets = 2, 3
>         for _ in range(history1_datasets):
>             self.dataset_populator.new_dataset(self.history1_id)
>         for _ in range(history2_datasets):
>             self.dataset_populator.new_dataset(self.history2_id)
> 
>         payload = { "history_id": self.history1_id, }
>         response = self._get("datasets", payload).json()
>         assert len(response) == history1_datasets
> 
>         payload = { "history_id": self.history2_id, }
>         response = self._get("datasets", payload).json()
>         assert len(response) == history2_datasets
> ```
{: .solution }
