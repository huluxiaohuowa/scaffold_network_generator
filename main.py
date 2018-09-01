# from multiprocessing import Pool,Manager
# from os import path
# import data
# import linecache

# num_process = 30
# input_file = 'input.txt'
# output_file = 'output.txt'
# input_dir = path.join(path.dirname(__file__),
#                       'data',
#                       'datasets',
#                       input_file
#                       )
#
# num_lines = -1
# for num_lines, line in enumerate(open(input_dir, 'r')):
#     pass
# num_lines += 1
#
#
# def get_smiles(idx_line):
#     return linecache.getline(input_dir, idx_line).strip()
#
#
# if __name__ == '__main__':
#     print(num_lines)
#     for i in range(num_lines):
#         try:
#             data.output_to_file(get_smiles(i+1))
#         except:
#             print(get_smiles(i+1))
# from data import *
#
# dic_scaffold, pro_scaffold = get_sng_protobuf()
import data
from multiprocessing import Manager, Pool
from os import path

manager = Manager()
q = manager.Queue()

p = Pool(processes=30)
for i in range(data.get_num_lines(path.join(path.dirname(__file__),
                                            'data',
                                            'datasets',
                                            'input.txt'
                                            ))):
    p.apply_async(data.sng_from_line_2_queue,
                  (i, q,)
                  )
p.close()
p.join()

scaffold_dict, dataset = data.protobuf_from_queue(q)
