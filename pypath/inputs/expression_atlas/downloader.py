
import aiohttp
import asyncio
import requests
import re
import tempfile
import zipfile
import os

import pypath.share.cache as cache
import pypath.share.session as session



from tqdm import tqdm



class ExpressionAtlasExperimentDownloader:
    def __init__(self) -> None:
        self.cache_dir = cache.get_cachedir()
        self._logger = session.Logger(name='inputs.expression_atlas.downloader')
        self._log = self._logger._log

    def download_files(self, urls: list[str], experiment_names: list[str], cache: bool = True) -> bool:
        """
        Downloads files from the given URLs and saves them with the corresponding experiment names.
        Args:
            urls (list[str]): A list of URLs to download files from Expression Atlas.
            experiment_names (list[str]): A list of experiment names corresponding to each experiement in Expression Atlas.
            See: https://www.ebi.ac.uk/gxa/experiments
        Returns:
            bool: True if all files are downloaded and saved successfully.
        Raises:
            ValueError: If the lengths of `urls` and `experiment_names` do not match.
            ValueError: If the file type is unsupported.
        Notes:
            - Handles ZIP files by calling the `_handle_zip_file` method.
            - Handles TSV files by calling the `_save_tsv_file` method.
            - Prints an error message if the download fails with a non-200 HTTP status code.
        """

        if len(urls) != len(experiment_names):
            raise ValueError("URLs and experiment names must have the same length")
        
        for url, experiment_name in tqdm(zip(urls, experiment_names), total=len(urls), desc="Downloading files"):
            response = requests.head(url)
            d = response.headers.get('content-disposition')
            fname = re.findall(r'filename="?([^"]+)"?', d)[0] if d else f"{experiment_name}.tsv"
            fname = fname.strip('\"').strip()
            
            if cache:
                self._handle_cache(fname)
                if self.use_cache:
                    continue
            
            response = requests.get(url)
            if response.status_code == 200:
                content_type = response.headers.get('Content-Type', '')
                
                # Handle ZIP file
                if "zip" in content_type or fname.endswith(".zip"):
                    self._handle_zip_file(response.content, experiment_name)
                # Handle TSV file
                elif "text/tab-separated-values" in content_type or fname.endswith(".tsv"):
                    self._save_tsv_file(response.content, fname)
                else:
                    raise ValueError(f"Unsupported file type {fname}")
            else:
                self._log(f"Download failed for {experiment_name}, HTTP Status Code: {response.status_code}")
        
        return True

    def _handle_cache(self, fname: str) -> None:
        """Handle caching for a given filename."""
        self.cache_file_name = os.path.join(self.cache_dir, fname)
        self.use_cache = os.path.exists(self.cache_file_name) and os.stat(self.cache_file_name).st_size > 0
        if self.use_cache:
            self._log('Cache file found, no need for download.')

    def _handle_zip_file(self, content: bytes, experiment_name: str) -> None:
        """
        Handles a ZIP file by extracting TSV files to a cache directory.
        Args:
            content (bytes): The content of the ZIP file in bytes.
            experiment_name (str): The name of the experiment, used to name the temporary ZIP file.
        Raises:
            FileNotFoundError: If no TSV file is found in the ZIP archive.
        """

        with tempfile.TemporaryDirectory() as tmpdirname:
            zip_file_path = os.path.join(tmpdirname, f"{experiment_name}.zip")
            with open(zip_file_path, 'wb') as zip_file:
                zip_file.write(content)
            
            with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
                found = False
                for file in zip_ref.namelist():
                    if file.endswith(".tsv"):
                        zip_ref.extract(file, self.cache_dir)
                        found = True
                        break
                
                if not found:                
                    raise FileNotFoundError(f"No TSV file found in the ZIP archive. Experiement: {experiment_name}")
        
    
    def _save_tsv_file(self, content: bytes, fname: str) -> None:
        """Save a TSV file to the cache directory."""
        tsv_file_path = os.path.join(self.cache_dir, fname)
        with open(tsv_file_path, 'wb') as tsv_file:
            tsv_file.write(content)



class AsyncExpressionAtlasExperimentDownloader(ExpressionAtlasExperimentDownloader):
    def __init__(self) -> None:
        super().__init__()
    
    async def download_files(self, urls: list[str], experiment_names: list[str]) -> bool:
        """
        Asynchronously downloads files from the given URLs and saves them with the corresponding experiment names.
        Args:
            urls (list[str]): A list of URLs to download files from Expression Atlas.
            experiment_names (list[str]): A list of experiment names corresponding to each experiement in Expression Atlas.
            See: https://www.ebi.ac.uk/gxa/experiments
        Returns:
            bool: True if all files are downloaded and saved successfully.
        Raises:
            ValueError: If the lengths of `urls` and `experiment_names` do not match.
            ValueError: If the file type is unsupported.
        Notes:
            - Handles ZIP files by calling the `_handle_zip_file` method.
            - Handles TSV files by calling the `_save_tsv_file` method.
            - Prints an error message if the download fails with a non-200 HTTP status code.
        """

        if len(urls) != len(experiment_names):
            raise ValueError("URLs and experiment names must have the same length")
        
        async with aiohttp.ClientSession() as session:
            tasks = []
            for url, experiment_name in zip(urls, experiment_names):
                tasks.append(self._download_file(session, url, experiment_name))
            
            for f in tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Downloading files"):
                await f
        
        return True

    async def _download_file(self, session: aiohttp.ClientSession, url: str, experiment_name: str, cache: bool = True) -> None:
        """
        Asynchronously downloads a file from the given URL and processes it based on its content type.

        Args:
            session (aiohttp.ClientSession): The aiohttp client session to use for making the request.
            url (str): The URL from which to download the file.
            experiment_name (str): The name of the experiment, used to generate the filename if not provided in the response.

        Raises:
            ValueError: If the file type is unsupported.

        Notes:
            - If the file is a ZIP file, it will be handled by the `_handle_zip_file` method.
            - If the file is a TSV file, it will be saved using the `_save_tsv_file` method.
            - If the HTTP status code is not 200, it will print an error message.
        """
        try:
            async with session.head(url) as response:
                d = response.headers.get('content-disposition')
                fname = re.findall(r'filename="?([^"]+)"?', d)[0] if d else f"{experiment_name}.tsv"
                fname = fname.strip('\"').strip()
            
            if cache:
                self._handle_cache(fname)
                if self.use_cache:
                    return

            async with session.get(url) as response:
                if response.status == 200:
                    content_type = response.headers.get('Content-Type', '')

                    # Handle ZIP file
                    if "zip" in content_type or fname.endswith(".zip"):
                        content = await response.read()
                        self._handle_zip_file(content, experiment_name)
                    # Handle TSV file
                    elif "text/tab-separated-values" in content_type or fname.endswith(".tsv"):
                        content = await response.read()
                        self._save_tsv_file(content, fname)
                    else:
                        raise ValueError(f"Unsupported file type {fname}")
                else:
                    self._log(f"Download failed, HTTP Status Code: {response.status}")

        except Exception as e:

            print(f"Async download error for {url}: {e}")

            self._log(
                f"Async download error for {url}: {e}"
            )
