from multiprocessing import Manager, cpu_count, Pool, Process
from data_utils import *
import argparse
from os import path


parser = argparse.ArgumentParser()
parser.add_argument(
    "--file_input",
    help="The location of input file",
    type=str,
    default=path.join(path.dirname(__file__),
            'data',
            'datasets',
            'input.txt')
    )
parser.add_argument(
    "--scaffolds_output",
    help="The location of scaffold SMILES file",
    type=str,
    default=path.join(path.dirname(__file__),
        'data',
        'datasets',
        'scaffolds.smi')
    )
parser.add_argument(
    "--file_output",
    help="The location of output file",
    type=str,
    default=path.join(path.dirname(__file__),
        'data',
        'datasets',
        'scaffolds.bin')
    )
parser.add_argument(
    "--np", 
    help="Num of processes", 
    type=int,
    default=cpu_count() - 1 
)
parser.add_argument(
    "--nq", 
    help="Num of items in the queue", 
    type=int
)
parser.add_argument(
    "--print_step", 
    help="Num of items in the queue", 
    type=int,
    default=5000
)
parser.add_argument(
    "--out_fmt",
    help="Format of output file",
    type=str,
    default="protobuf"
)
parser.add_argument(
    "--sql_dic",
    help="A json file with postgres connection kargs",
    type=str,
    default='data1.json'
)
parser.add_argument(
    "--db_name",
    help="Database name",
    type=str,
    default='sc1'
)

args = parser.parse_args()
q = Manager().Queue(args.nq)
q.put(args.file_input)
p = Pool(processes=args.np)

if args.out_fmt == "protobuf":
    p_get = Process(
        target=data_from_queue,
        args=(
            q,
            args.file_output,
            args.scaffolds_output,
            args.print_step
        )
    )
elif args.out_fmt == "postgres":
    p_get = Process(
        target=sql_from_queue,
        args=(
            q,
            args.sql_dic,
            args.db_name,
            args.print_step
        )
    )
else:
    raise ValueError("out_fmt must be \"protobuf\" or \"postgres\"")


for i in range(get_num_lines(args.file_input)):
    p.apply_async(
        sng_from_line_2_queue,
        (
            i,
            q,
            args.file_input
        )
    )


p_get.start()
p.close()
p.join()
p_get.join()
p_get.close()
# sng_to_queue(q, processes=args.np, file=args.file_input)



# dic_scaffold, dic_idx_sm = data_from_queue(q)

# with open(args.file_output, 'wb') as f:
#     f.write(dic_scaffold.SerializeToString())
# with open(args.scaffolds_output, 'wb') as f:
#     f.write(dic_idx_sm.SerializeToString())
