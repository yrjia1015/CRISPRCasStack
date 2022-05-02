import os

def delete_folderfile(path):
    file_list = os.listdir(path)
    for file in file_list:
        # print(file)
        file_path=os.path.join(path,file)
        os.remove(file_path)
    return 0