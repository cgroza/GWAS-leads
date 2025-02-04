import pandas as pd
import gzip
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Get lead Commands:')
    parser.add_argument('--raw', type=str, default=None, help='path to summary statistics', required=True)
    parser.add_argument('--save', type=str, default=None, help='path to save results', required=True)
    parser.add_argument('--gwp', type=float, default=5e-8, help='genome-wide significant p-value')
    parser.add_argument('--dist', type=int, default=1e6, help='distance between leads')
    parser.add_argument('--pcol', type=str, default=None, help='p-value', required=True)
    parser.add_argument('--chr', type=str, default=None, help='CHR', required=True)
    parser.add_argument('--pos', type=str, default=None, help='POS', required=True)
    parser.add_argument('--a1', type=str, default=None, help='A1', required=False)
    parser.add_argument('--a2', type=str, default=None, help='A2', required=False)
    args = parser.parse_args()
    return args

def get_lead(args):
    leads = []
    target = None
    distance = int(args.dist)
    with gzip.open(args.raw, 'rt') as ss:
        header = next(ss)
        for line in ss:
            ll = dict(zip(header.strip().split(), line.strip().split()))
            # filter p-value, keep SNP, remove HLA
            if float(ll[args.pcol]) < args.gwp:   # and len(ll[args.a1])==len(ll[args.a2])==1:
                # if ll[args.chr]=='6' and int(ll[args.pos])<34448354 and int(ll[args.
                #     continue
                if target is None:
                    target = ll
                else:
                    if ll[args.chr]!=target[args.chr] or int(ll[args.pos]) - int(target[args.pos]) > distance:
                        leads.append(target)
                        target = ll
                    else:
                        if float(ll[args.pcol]) < float(target[args.pcol]):
                            target = ll
    leads.append(target)
    return leads

if __name__ == '__main__':
    args = parse_args()
    leads = get_lead(args)
    df_lead = pd.DataFrame(leads)
    df_lead.to_csv(args.save, sep='\t', index=False, header=True)
    df_gtf = df_lead[[args.chr]].copy()
    df_gtf.loc[:, 'start'] = (df_lead[args.pos].astype(int) - args.dist/2).astype(int)
    df_gtf.loc[:, 'end'] = (df_lead[args.pos].astype(int) + args.dist/2).astype(int)
    df_gtf.loc[df_gtf['start'] < 0, 'start'] = 0
    df_gtf.to_csv('{}.gtf'.format(args.save), sep='\t', index=False, header=False)
