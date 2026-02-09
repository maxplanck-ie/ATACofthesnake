import hashlib
from pathlib import Path
from rich.console import Console
from aos.logger import setup_logger
import requests
import click
from importlib.metadata import version

@click.command(
        context_settings=dict(
            help_option_names=["-h", "--help"]
            )
    )
@click.option(
    '-o',
    '--outputdir',
    default='aos_example',
    help='Specify a directory to download the example data to.'
)
@click.argument(
    'dataset',
    default='both',
    type=click.Choice(['both', 'immgen', 'cele'], case_sensitive=False),
)
def download_example(dataset: str, outputdir: str):
    AOSExample(outputdir, dataset)


class AOSExample:
    '''
    Download example data for AOS.
    '''

    def __init__(self, output: str, dataset: str) -> None:
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        self.dataset = dataset
        self.base_url = 'https://zenodo.org/record/18553480/files/'
        self.wanted_files = {
            'cele': {
                'cele_accessions.txt': 'e5dc7bdc2af68fd2d59b8ac6bae99bcd',
                'cele_genes.gtf': '38579d77ef01a295b7009dfcfdc2c32b',
                'cele_genome.fa': 'd570defcdc006a7c2859fc92dbb21bc4',
                'cele_genome.fa.fai': '024293aa8ce44cdaa5651fcf3d7492a5',
                'cele_rar.bed': '4d5046d8914a2d60757577df38daf92e',
                "Auxin_10h_rep1.cram": "9c1d54fca0eae86a183c33c0962557f9",
                "Auxin_10h_rep2.cram": "30fc30607c5d12098dea337e192cd7ea",
                "Auxin_11h_rep1.cram": "4c1cdd0c867cbd8be0472f14dfad5c46",
                "Auxin_1h_rep1.cram": "e43b5e7d361ae52fdf9b33db6202f94b",
                "Auxin_1h_rep2.cram": "3a21fa9c12f3cf995bd3c93ad4cc406c",
                "Auxin_2h_rep1.cram": "4e40293b6f7769cfe0f2f8906bd0d71b",
                "Auxin_2h_rep2.cram": "64416558b15fae7bbf85d1b3158e16bc",
                "Auxin_3h_rep1.cram": "13a0f89606c63e4d2c1bbdb4a0dc5ed0",
                "Auxin_3h_rep2.cram": "13679f1988e9e827b632b320c9bdfecb",
                "Auxin_4h_rep1.cram": "6eedcd2bdb500d7e261566cfc5cd247f",
                "Auxin_4h_rep2.cram": "321a7a5ecdf250be03efe7ae56a1c47f",
                "Auxin_5h_rep1.cram": "3afcd6c4e068e26e34ba0f024f5dbde6",
                "Auxin_5h_rep2.cram": "bfdc88678463d3a0d53df17fa74bfafb",
                "Auxin_6h_rep1.cram": "bf9c05bcf67db7632b0ed0cb1a2aaf89",
                "Auxin_6h_rep2.cram": "3668663d3d483198e091ebf47f41ba6a",
                "Auxin_7h_rep1.cram": "7623b876583cb2560a95b00342bfd1e4",
                "Auxin_7h_rep2.cram": "bf6d8401f111c49521a86d45269dc990",
                "Auxin_8h_rep1.cram": "6ba00cedb88143c8a056415cdfd82926",
                "Auxin_8h_rep2.cram": "61a3d0fb2d9b14c3f542069a9b3f1b19",
                "Auxin_9h_rep1.cram": "9e5d18cb32c6484c332355031adec043",
                "Auxin_9h_rep2.cram": "a5e0184cc19591f7579ff17ead83a027",
                "EtOH_10h_rep1.cram": "1d48e3682df73bafb7038262867c3ea0",
                "EtOH_10h_rep2.cram": "9578175c0a4d9bea6cfce73d8263bab5",
                "EtOH_11h_rep1.cram": "4b1877597b2fe4cc73f6626a96417448",
                "EtOH_1h_rep1.cram": "8adad133373aea86293cd7053074e4c3",
                "EtOH_1h_rep2.cram": "1697f24c41b110929ec2fed93a973952",
                "EtOH_2h_rep1.cram": "79b769eaeeb4698d218936f15fbd5fd8",
                "EtOH_2h_rep2.cram": "fea3492214d1985e57d7df6ea48b2f13",
                "EtOH_3h_rep1.cram": "3af1f10b0a5c1ff098454587e7939903",
                "EtOH_3h_rep2.cram": "c9cb3ddabefead7db57199781e57a490",
                "EtOH_4h_rep1.cram": "af2eb56f6291663da031bc2de9869667",
                "EtOH_4h_rep2.cram": "bc86a2ec804f45fa8282b0d8d48b752a",
                "EtOH_5h_rep1.cram": "27fef1210576a91c8a2bb790afefb251",
                "EtOH_5h_rep2.cram": "68b1a2238e51188d28068cc0bce91649",
                "EtOH_6h_rep1.cram": "c30333b6e1426f888f496a0b4c298abb",
                "EtOH_6h_rep2.cram": "231d7df24436cb1fc5a4efb30007963c",
                "EtOH_7h_rep1.cram": "2a5cafd47a1928d5ad08f9bb35c9eb98",
                "EtOH_7h_rep2.cram": "bbb35b8682d80efe3c2bab934d843205",
                "EtOH_8h_rep1.cram": "8ead8d3543acacca7516af4a8897176e",
                "EtOH_8h_rep2.cram": "d953c5148b74019e92f91b4f84a80b3f",
                "EtOH_9h_rep1.cram": "f68fd6da0f241446117460623432e85f",
                "EtOH_9h_rep2.cram": "6959a75f03f85e059e5584f7ee430a5e",
                "Untreated_0h_rep1.cram": "825b712fd8ff3c8ce0ba5e7fb11aba65"               
            },
            'immgen': {
                'immgen_accessions.txt': '64bd99e8fb4098041d8068319b1671ef',
                'immgen_genes.gtf': '1b63667effee9e06b46a189a126bc5f6',
                'immgen_genome.fa': '32cfe205a8fc667094334f74085dac80',
                'immgen_genome.fa.fai': 'e975a9f9949b3de65167b1bd82352952',
                'immgen_rar.bed': '0cd3de126237f3d6d9da9ef023caadf7',
                "preT_DN1_Th_rep1.cram": "8d833cbfb37351e74f9a69566436f245",
                "preT_DN1_Th_rep2.cram": "364a2e0b5ba629cd233b666bcf0fe084",
                "preT_DN2b_Th_rep1.cram": "d8db733101847b4dda7c142cd18c7eae",
                "preT_DN2b_Th_rep2.cram": "022aefa27a8ecf71f3663e86bc50c5b9",
                "preT_DN3_Th_rep1.cram": "38cbacb48a79b6dab9d67f382276e5b5",
                "preT_DN3_Th_rep2.cram": "0f0fe0043cb08d777bf9e99f129967fd",
                "T_4_Th_rep1.cram": "205a0731dc549554bfe3dcbd846e88b5",
                "T_4_Th_rep2.cram": "f00aefbda0b0839ad21e6eebd8b5b538",
                "T_8_Th_rep1.cram": "eed0f8f9e42c5bae95e26044a9a9d852",
                "T_8_Th_rep2.cram": "837003be71dd5614b5db899b7d039123",
                "T_DP_Th_rep1.cram": "537d90b9d86158501eeb8daffa8a511a",
                "T_DP_Th_rep2.cram": "194cc34b62b68e49f6842de6941f950c",
            }
        }


        # Set up log
        self.logfile = self.output / 'AOS_example_download.log'
        self.logger = setup_logger(self.logfile)
        self.logger.info(f"AOS - version {version('ATACofthesnake')}")
        self.download_and_verify()


    def download_and_verify(self):
        '''
        Download the actual data.
        '''
        self.logger.info(f"Downloading example data from zenodo to {self.output.resolve()}")
        match self.dataset:
            case 'both':
                wanted = [(dataset, fname, md5) for dataset, files in self.wanted_files.items() for fname, md5 in files.items()]
                Path(self.output / 'immgen').mkdir(exist_ok=True)
                Path(self.output / 'cele').mkdir(exist_ok=True)
            case 'immgen':
                wanted = [(dataset, fname, md5) for dataset, files in self.wanted_files.items() if dataset == 'immgen' for fname, md5 in files.items()]
                Path(self.output / 'immgen').mkdir(exist_ok=True)
            case 'cele':
                wanted = [(dataset, fname, md5) for dataset, files in self.wanted_files.items() if dataset == 'cele' for fname, md5 in files.items()]
                Path(self.output / 'cele').mkdir(exist_ok=True)
        
        verified = 0
        for file in wanted:
            dataset, fname, md5 = file
            url = self.base_url + fname
            if (self.output / dataset / fname).exists():
                self.logger.info(f"File {fname} already exists, verifying md5. File {verified + 1} of {len(wanted)}")
                md5_hash = hashlib.md5()
                with open(self.output / dataset / fname, "rb") as f:
                    # Read and update hash string value in blocks of 4K
                    for byte_block in iter(lambda: f.read(4096), b""):
                        md5_hash.update(byte_block)
                if md5_hash.hexdigest() != md5:
                    self.logger.warning(f"MD5 checksum does not match for existing file {fname}, re-downloading.")
                else:
                    self.logger.info(f"MD5 checksum matches for existing file {fname}, skipping download.")
                    verified += 1
                    continue
            self.logger.info(f"Downloading {fname} from {url}...")
            response = requests.get(url)
            with open(self.output / dataset / fname, 'wb') as f:
                f.write(response.content)
            self.logger.info(f"Downloaded {fname}, verifying md5. File {verified + 1} of {len(wanted)}")
            md5_hash = hashlib.md5()
            with open(self.output / dataset / fname, "rb") as f:
                # Read and update hash string value in blocks of 4K
                for byte_block in iter(lambda: f.read(4096), b""):
                    md5_hash.update(byte_block)
            if md5_hash.hexdigest() != md5:
                self.logger.error(f"MD5 checksum does not match for {fname}, Exiting.")
                raise ValueError(f"MD5 checksum does not match for {fname}, something went wrong during download.")
            self.logger.info(f"MD5 checksum matches for {fname}.")
            verified += 1
        self.logger.info(f"All files downloaded and verified successfully under {self.output}")