import subprocess
from cookiecutter.main import cookiecutter

def main():
    cookiecutter(
        'https://github.com/yztxwd/snakemake-pipeline-general.git',
        extra_context={'directory_name': 'snakemake-pipeline-general'}
    )
    
if __name__ == '__main__':
    main()
