> ### {% icon comment %} This is best viewed in Pangeo Jupyter Lab
>
> 1. Create a history
>
>     > ### {% icon hands_on %} Hands-on: Create history
>     >
>     > 1. Make sure you start from an empty analysis history.
>     > 2. **Rename your history** to be meaningful and easy to find. For instance, you can choose **Xarray with Pangeo notebook** as the name of your new history.
>     {: .hands_on}
>
>
> 2. Upload CAMS PM2.5 data
>
>    Data can be retrieved directly from [Copernicus Atmosphere Monitoring Service](https://ads.atmosphere.copernicus.eu/) but to make it easier, you can download the tutorial data from [Zenodo](https://zenodo.org/record/5805953/files/CAMS-PM2_5-20211222.netcdf).
>
>    > ### {% icon hands_on %} Hands-on: Data upload
>    >
>    > 1. Create a new history for this tutorial
>    > 2. Import the files from [Zenodo](https://doi.org/10.5281/zenodo.5805953) or from
>    >    the shared data library (`GTN - Material` -> `Climate` -> `Pangeo Notebook in Galaxy - Introduction to Xarray`):
>    >
>    >    ```
>    >    https://zenodo.org/record/5805953/files/CAMS-PM2_5-20211222.netcdf
>    >    ```
>    >
>    > 3. Rename the dataset to **CAMS-PM2_5-20211222.netcdf**
>    > 4. Check that the datatype is **netcdf**
>    >
>    {: .hands_on}
>
> 3. Starting Galaxy Pangeo JupyterLab
>
>    > ### {% icon hands_on %} Hands-on: Launch Pangeo notebook JupyterLab in Galaxy
>    >
>    > Currently Pangeo notebook in Galaxy is available on [useGalaxy.eu](https://usegalaxy.eu) only. Make sure you login to Galaxy and search for Pangeo notebook and not the default JupyterLab in Galaxy to make sure you ahve all the Pangeo Software stack available. The default JupyterLab in Galaxy would not be sufficient for executing all the tasks in this tutorial.
>    >
>    > 1. Open the {% tool [Pangeo Notebook](interactive_tool_pangeo_notebook) %} by clicking [here](https://usegalaxy.eu/?tool_id=interactive_tool_pangeo_notebook){:target="_blank"}
>    > 2. *Include data into the environment*: `CAMS-PM2_5-20211222.netcdf`
>    > 3. Click on **Execute**: the tool will start running and will stay running
>    > 4. Click on the **User** menu at the top and go to **Active Interactive Tools** and locate the Pangeo JupyterLab instance you started.
>    > 5. Click on your Pangeo JupyterLab instance
>    {: .hands_on}
>
> 4. Once the notebook has started, open a Terminal in Jupyter **File → New → Terminal**
>
> 5. Run the following command:
>
>    ```
>    wget {{ ipynbpath }}
>    ```
>
> 6. Switch to the Notebook that appears in the list of files on the left
{: .comment}
