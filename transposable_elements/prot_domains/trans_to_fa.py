import sys

def process_file(filename, genomeid):
    with open(filename, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            line=line.strip().split('\t')
            if line_number % 2 ==1:
                print(">"+genomeid+'_'+line[0]+'-'+line[5]+':'+line[6]+'..'+line[7])
            else:
                print(line[1])


if __name__ == "__main__":
    filename=sys.argv[1]
    genomeid=sys.argv[2]
    process_file(filename, genomeid)
