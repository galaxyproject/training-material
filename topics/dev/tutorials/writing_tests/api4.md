
> <solution-title></solution-title>
> 
> Using populators.
> 
> ```python
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
