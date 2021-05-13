
Changes to file ``lib/galaxy_test/api/test_users.py``:

```diff
index ca1796e..47b24e6 100644
--- a/lib/galaxy_test/api/test_users.py
+++ b/lib/galaxy_test/api/test_users.py
@@ -171,6 +171,33 @@ class UsersApiTestCase(ApiTestCase):
         search_response = get(url).json()
         assert "cat1" in search_response
 
+    def test_favorite_extensions(self):
+        index_response = self._get("users/current/favorites/extensions")
+        index_response.raise_for_status()
+        index = index_response.json()
+        assert isinstance(index, list)
+        assert len(index) == 0
+
+        create_response = self._post("users/current/favorites/extensions/fasta")
+        create_response.raise_for_status()
+
+        index_response = self._get("users/current/favorites/extensions")
+        index_response.raise_for_status()
+        index = index_response.json()
+        assert isinstance(index, list)
+        assert len(index) == 1
+
+        assert "fasta" in index
+
+        delete_response = self._delete("users/current/favorites/extensions/fasta")
+        delete_response.raise_for_status()
+
+        index_response = self._get("users/current/favorites/extensions")
+        index_response.raise_for_status()
+        index = index_response.json()
+        assert isinstance(index, list)
+        assert len(index) == 0
+
     def __url(self, action, user):
         return self._api_url("users/{}/{}".format(user["id"], action), params=dict(key=self.master_api_key))
 
-- 
2.7.4


```
