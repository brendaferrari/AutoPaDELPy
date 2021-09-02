import pandas as pd
import os
import shutil
import glob

def moveOutputFile():
    print('Moving the output files to output_files directory...')
    source_folder = os.getcwd()
    directory = 'output_files'
    path = os.path.join(source_folder, directory)

    os.mkdir(path)
    destination_folder = path

    # to get files to move
    files_to_move = ('*.csv', '*.out')
    files_grabbed = []
    for files in files_to_move:
        files_grabbed.extend(glob.glob(files))

    # iterate files
    for file in files_grabbed:
        # construct full file path
        source = source_folder + '/' + file # Adding *, try to take the * out
        destination = destination_folder + '/' + file
        # move file
        shutil.move(source, destination)
        print('Moved:', file)
    return 

def saveFile():

    path = os.getcwd() + '/' + 'output_files'
    all_files = glob.glob(path + "/*.csv")

    li = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        li.append(df)

    frame = pd.concat(li, axis=1)
    frame.to_csv('output_files/descriptors_unprocessed.csv')

    print(f'The files were saved inside output_files directory!')
    return