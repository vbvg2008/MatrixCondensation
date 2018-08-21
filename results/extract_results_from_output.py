import numpy as np
import os


def read_and_decode_mcge_p(path,file):
    with open(os.path.join(path,file), 'r') as f:
        content = f.read()
    number = []
    for s in content.split():
        try:
            x = float(s)
            number.append(x)
        except:
            pass
    content_matrix = np.reshape(number, (8, 5, 8))
    num_proc = content_matrix[0, 0, 2]
    average_time = content_matrix[:,:,3]
    comm_time = content_matrix[:,:,5]
    dist_time = content_matrix[:,:,7]
    return num_proc, average_time, comm_time, dist_time

def read_and_decode_ge_scalapack(path, file):
    with open(os.path.join(path,file), 'r') as f:
        content = f.read()
    content_list = content.split()
    content_2 = content_list[9::9]
    number = []
    for content_2_element in content_2:
        for s in content_2_element.split(','):
            try:
                x = float(s)
                number.append(x)
            except:
                pass
    content_matrix = np.reshape(number, (8, 5, 6))
    average_time = content_matrix[:, :, 5]
    tmp_number = content_list[8].split(',')
    num_proc = int(tmp_number[0])
    return num_proc, average_time



def generate_metrics_mcge_p(path):
    '''
    this function will generate 3x8x8x5 matrices for mc_p and ge_p
    1st dimension: 0-total execution time  1-total communication time 2-data distribution time
    2nd dimension:problem size in thousand (1-8K)
    3rd dimension: number of processors-1 2 4 8 16 32 64 128
    4th dimension: different independent runs'''
    
    #initialize result
    result = np.zeros((3,8,8,5), dtype = float)
    file_list = os.listdir(path)
    for file in file_list:
        num_proc, average_time, comm_time, dist_time = read_and_decode_mcge_p(path,file)
        index = int(np.log2(num_proc))
        result[0, :, index, :] = average_time
        result[1, :, index, :] = comm_time
        result[2, :, index, :] = dist_time
    return result

def generate_metrics_gescalapack(path):
    '''
    This function will generate 8x8x5 matrices for ge_scalapack
    '''
    #initialize result
    result = np.zeros((8, 8, 5), dtype = float)
    file_list = os.listdir(path)
    for file in file_list:
        num_proc, average_time = read_and_decode_ge_scalapack(path,file)
        index = int(np.log2(num_proc))
        result[:, index, :] = average_time
    return result


if __name__ == "__main__":
    
    mc_dir = "../output/mc_p"
    ge_dir = "../output/ge_p"
    gescalapack_dir = "../output/ge_scalapack"
    result_mc = generate_metrics_mcge_p(mc_dir)
    result_ge = generate_metrics_mcge_p(ge_dir)
    result_gescalapack = generate_metrics_gescalapack(gescalapack_dir)
    
    #averaging 5 different runs
    result_mc_average = np.mean(result_mc, axis = -1)
    result_gep_average = np.mean(result_ge, axis = -1)
    result_gescalapack_average = np.mean(result_gescalapack, axis = -1)

