
> <solution-title></solution-title>
> 
> API test for creating an object.
> 
> ```python
> class MyTutorialApiTestCase(ApiTestCase):
> 
>     def test_create_role(self):
>         # prepare new role
>         name = "cool role"
>         description = "description of this cool role"
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
> ```
{: .solution }
