import os
import shutil
import subprocess

CHILD_DIR = 'snakemake-pipeline-general/'

def main():
    # concatenate config file
    with open(f'{CHILD_DIR}/config.yaml') as fi, open('config.yaml', 'a') as fo:
        shutil.copyfileobj(fi, fo)
    
    # include rules into Snakefile
    with open('Snakefile', 'a') as f:
        for rule in os.listdir(f'{CHILD_DIR}/rules'):
            f.write(f'include: {CHILD_DIR}/rules/{rule}\n')
        for rule in os.listdir('./rules'):
            f.write(f'include: ./rules/{rule}\n')
    
    # init the git
    subprocess.check_call(['git', 'init'], stdout=subprocess.DEVNULL)
    subprocess.check_call(['git', 'add', '.'], stdout=subprocess.DEVNULL)
    subprocess.check_call(['git', 'commit', '-m', '"first commit"'], stdout=subprocess.DEVNULL)

if __name__ == '__main__':
    main()
