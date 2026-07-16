import pynmrstar
import os
from datetime import date
import plotly.express as px
import logging
logging.getLogger('pynmrstar').setLevel(logging.ERROR)


from matplotlib.pyplot import xlabel, ylabel


def remediate_relaxation_data(str_file):
    pass

def check_entry_for_relaxation_data(str_file):
    ent = pynmrstar.Entry.from_file(str_file,raise_parse_warnings=False)
    seq = ent.get_tag('Entity.Polymer_seq_one_letter_code')[0].replace("\n", "")
    print(ent.entry_id, seq)
    print(calculate_molecular_weight(seq))
    remediated_ent = pynmrstar.Entry.from_scratch(ent.entry_id)
    author_provided_unit_t1=[]
    author_provided_unit_t2=[]
    s2=0
    t1rho=0
    hetnoe=0
    field_strength=[]
    for saveframe in ent:
        # if saveframe.category == 'heteronucl_T1_relaxation':
        #     t1+=1
        # if saveframe.category == 'heteronucl_T2_relaxation':
        #     t2+=1
        if saveframe.category =='heteronucl_NOEs':
            hetnoe+=1
            field_strength.append(saveframe.get_tag('Spectrometer_frequency_1H')[0])
        if saveframe.category == 'order_parameters':
            s2+=1
        if saveframe.category == 'heteronucl_T1rho_relaxation':
            t1rho+=1
        if saveframe.category in ['heteronucl_T1_relaxation','heteronucl_T2_relaxation']:
            # units found in entries [['s'], ['ms'], ['s-1'], ['ms-1'], ['ns'], ['Hz']]
            if saveframe.category == 'heteronucl_T1_relaxation':
                author_provided_unit = saveframe.get_tag('T1_val_units')[0]
                author_provided_unit_t1.append(saveframe.get_tag('T1_val_units')[0])
                saveframe.add_tag(name='T1_val_units',value='s', update=True)
                if ent.entry_id == '51478':
                    author_provided_unit = 's-1' #verified in the manuscript
                if ent.entry_id == '17701':
                    author_provided_unit = 'ms' # based on the manuscript
                if ent.entry_id in ['51305','51306','51307']:
                    author_provided_unit = 's-1'
                if ent.entry_id == '16737':
                    author_provided_unit = 'ms' # wrongly mentioned as ms-1
                if ent.entry_id == '27646':
                    author_provided_unit = 's' # worngly mentioned as s-1
                if ent.entry_id in ['6881','51126','18231','18758','18087','6880','51224','51223']:
                    author_provided_unit = 's-1' # worngly mentioned as s
                if ent.entry_id in ['5154','7088','7056','25525','25523','30834']:
                    author_provided_unit = 's'
                saveframe.add_tag(name='T1_val_units',value='s', update=True)
            else:
                author_provided_unit = saveframe.get_tag('T2_val_units')[0]
                if ent.entry_id == '50734':
                    author_provided_unit = 'ms' # There may be type in T2 unit in the file
                if ent.entry_id == '51478':
                    author_provided_unit = 's-1'#verified in the mauscript
                if ent.entry_id == '17701':
                    author_provided_unit = 'ms'
                if ent.entry_id == '50745':  # wrongly mentioned as s
                    author_provided_unit = 'ms'
                if ent.entry_id in ['51305','51306','51307','18758']:
                    author_provided_unit = 's-1'
                if ent.entry_id == '16737':
                    author_provided_unit = 'ms' # wrongly mentioned as ms-1
                if ent.entry_id == '27646':
                    author_provided_unit = 's' # worngly mentioned as s-1
                if ent.entry_id in ['6881','51126','18231','18087','6880','51224','51223']:
                    author_provided_unit = 's-1' # worngly mentioned as s
                if ent.entry_id in ['5154','7088','7056','25525','25523','30834']:
                    author_provided_unit = 's'
                author_provided_unit_t2.append(saveframe.get_tag('T2_val_units')[0])
                saveframe.add_tag(name='T2_val_units',value='s', update=True)
            for loop in saveframe:
                if loop.category in ['_T1','_T2']:
                    for row in loop.data:
                        if loop.category == '_T1':

                            if author_provided_unit == 's':
                                pass
                            elif author_provided_unit == 'ms':

                                try:
                                    row[loop.tags.index('Val_err')] = round(float(row[loop.tags.index('Val_err')])/1000.00,3)
                                except ValueError:
                                    pass
                                row[loop.tags.index('Val')] = round(float(row[loop.tags.index('Val')]) / 1000.00, 3)
                            elif author_provided_unit in ['s-1','Hz']:
                                try:
                                    row[loop.tags.index('Val_err')] = round(
                                        float(row[loop.tags.index('Val_err')]) / pow(float(row[loop.tags.index('Val')]), 2), 3)
                                    try:
                                        row[loop.tags.index('Val')] = round(1.0 / float(row[loop.tags.index('Val')]), 3)
                                    except ZeroDivisionError:
                                        pass
                                except ValueError:
                                    pass
                                except ZeroDivisionError:
                                    pass
                            elif author_provided_unit == 'ms-1':

                                row[loop.tags.index('Val_err')] = round((float(row[loop.tags.index('Val_err')])/pow(float(row[loop.tags.index('Val')]),2))/1000.0 , 3)
                                row[loop.tags.index('Val')] = round(1000.0 / float(row[loop.tags.index('Val')]), 3)
                            elif author_provided_unit == 'ns':
                                try:
                                    row[loop.tags.index('Val_err')] = round(float(row[loop.tags.index('Val_err')])/1000000000.00,3)
                                    row[loop.tags.index('Val')] = round(
                                        float(row[loop.tags.index('Val')]) / 1000000000.00, 3)
                                except ValueError:
                                    pass
                            else:
                                raise ValueError(f'New unit found{author_provided_unit}')
                        else:
                            if author_provided_unit == 's':
                                pass
                            elif author_provided_unit == 'ms':
                                try:
                                    row[loop.tags.index('T2_val_err')] = round(float(row[loop.tags.index('T2_val_err')])/1000.00,3)
                                except ValueError:
                                    pass
                                row[loop.tags.index('T2_val')] = round(float(row[loop.tags.index('T2_val')]) / 1000.00,
                                                                       3)
                            elif author_provided_unit in ['s-1','Hz']:
                                try:

                                    row[loop.tags.index('T2_val_err')] = round(float(row[loop.tags.index('T2_val_err')])/pow(float(row[loop.tags.index('T2_val')]),2) , 3)
                                    row[loop.tags.index('T2_val')] = round(1.0 / float(row[loop.tags.index('T2_val')]),
                                                                           3)
                                except ValueError:
                                    pass
                                except ZeroDivisionError:
                                    pass
                            elif author_provided_unit == 'ms-1':
                                row[loop.tags.index('T2_val_err')] = round((float(row[loop.tags.index('T2_val_err')])/pow(float(row[loop.tags.index('T2_val')]),2))/1000.0 , 3)
                                row[loop.tags.index('T2_val')] = round(1000.0 / float(row[loop.tags.index('T2_val')]),
                                                                       3)
                            elif author_provided_unit == 'ns':
                                row[loop.tags.index('T2_val_err')] = round(float(row[loop.tags.index('T2_val_err')])/1000000000.00,3)
                                row[loop.tags.index('Val')] = round(float(row[loop.tags.index('Val')]) / 1000000000.00,
                                                                    3)
                            else:
                                raise ValueError(f'New unit found{author_provided_unit}')
        remediated_ent.add_saveframe(saveframe)
    #if len(author_provided_unit_t1)==0:
    if len(author_provided_unit_t1)!=0 and len(author_provided_unit_t2)!=0:
        print (ent.entry_id,len(author_provided_unit_t1)+len(author_provided_unit_t2),len(author_provided_unit_t1),len(author_provided_unit_t2))
        #print (ent.entry_id,len(author_provided_unit_t1),set(author_provided_unit_t1),len(author_provided_unit_t2),set(author_provided_unit_t2),hetnoe,s2,t1rho)

    #print (ent.entry_id,len(author_provided_unit_t1),author_provided_unit_t1,len(author_provided_unit_t2),author_provided_unit_t2,hetnoe,field_strength,t1rho,s2)
    #if set(author_provided_unit_t1) != set(author_provided_unit_t2):
      #  print (ent.entry_id,len(author_provided_unit_t1),set(author_provided_unit_t1),len(author_provided_unit_t2),set(author_provided_unit_t2))
    release_loop=remediated_ent.get_loops_by_category('_Release')
    #release_loop.add_data()
    comment= ''
    if ent.entry_id == '17701':
        comment = 'values were reported in ms but unit mentioned as s'
    if ent.entry_id == '19356':
        comment = 'T1 was okay, but T2 was given ms'
    if ent.entry_id == '50553':
        comment = 'poor data, some values are negative, manuscript not published, preprint available'
    if ent.entry_id == '27594':
        comment = 'one of the T1 value is negative, values are not reproted in the manuscript'
    if ent.entry_id == '50745':
        comment = 'T2 values were wrongly mentioned as s'
    if ent.entry_id in ['51305','51306','51307']:
        comment = 'author reported rates with with wrong units (T1 ms, T2 s)'
    if ent.entry_id == '16737':
        comment = 'wrongly mentioned as ms-1'
    if ent.entry_id == '27646':
        comment = 'values looks okay but looks like unit was mentioned worngly as s-1'
    if ent.entry_id in ['6881','51126','18231','18758','18087','6880']:
        comment = 'reported rates in the manuscript, but provided s as unit'
    if ent.entry_id in ['51224','51223']:
        comment = 'reported rates in the manuscript, but provided ms as unit'
    if ent.entry_id in ['5154','7088','7056','25525','25523']:
        comment = 'reported times , but provided s-1 as unit'
    if ent.entry_id in ['30834']:
        comment = 'reported in s , but provided ms as unit'
    release_loop[0].add_data({
        'Release_number':[str(max([int(i) for i in release_loop[0].get_tag('Release_number')])+1)],
        'Submission_date':[release_loop[0].get_tag('Submission_date')[0]],
        'Release.Type':['update'],
        'Author':['BMRB'],
        'Detail':[f'Relaxation (T1 & T2) data converted into s;{comment}'],
        'Entry_ID':[release_loop[0].get_tag('Entry_ID')[0]],
                             'Date':[str(date.today())]
    })
    remediated_ent.normalize()
    remediated_ent.write_to_file(f'/Users/kumaranbaskaran/Projects/bmrb/BMRB-data-remediation/Relaxation_data/data/output/bmr{ent.entry_id}_3.str')


def plot_data(directory_path):
    data={}
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)

        if os.path.isfile(file_path):  # Check if it's a file, not a directory
            data.update(extract_data(file_path))
            print (len(data))
    t1=[]
    dt1=[]
    t2=[]
    dt2=[]
    t=[]
    flg=[]
    good=[]
    bad=[]
    mw=[]
    unit=[]
    for k in data:
        if 'T1' in data[k] and 'T2' in data[k]:
            if data[k]['T1'] < data[k]['T2']:
                c = 'T1<=T2'
                if k[0] not in bad:
                    bad.append(k[0])
            else:
                c = 'T2>T2'
                if k[0] not in good:
                    good.append(k[0])
            unit.append(f'{data[k]['t1unit']},{data[k]['t2unit']}')
            t1.append(data[k]['T1'])
            mw.append(data[k]['mw'])
            dt1.append(data[k]['T1_err'])
            t2.append(data[k]['T2'])
            dt2.append(data[k]['T2_err'])
            t.append(k)
            flg.append(c)
    c1 = flg.count('T1<=T2')
    c2 = len(flg) - c1
    fig = px.scatter(x=t1, y=t2, hover_name=t, labels={'x': 'T1', 'y': 'T2','color':'T1 unit ,T2 unit'}, color=unit)
    fig.update_layout(
        xaxis=dict(range=[-100, 3000]),  # Set the x-axis range from 0 to 10
        # yaxis=dict(range=[-5, 5])  # Set the y-axis range from -5 to 5
    )
    fig.show()
    fig.write_html('T1_T2_units.html')
    fig.write_image('T1_T2_units.png', width=1200, height=800, scale=2)
    fig.write_image('T1_T2_units.jpg', width=1200, height=800)
    fig.write_image('T1_T2_units.pdf', width=1200, height=800)
    fig = px.scatter(x=t1,y=t2,hover_name=t,labels = {'x':'T1[?]', 'y' : 'T2[?]'},color=flg,title=f'T1<T2={c1}({len(bad)}) T1>T2={c2}({len(good)})')
    fig.write_html('before_correction.html')
    fig.show()
    fig2 = px.scatter(x=mw, y=t1, hover_name=t, labels={'x': 'Molecular Weight', 'y': 'T1[s]'}, color=flg,
                      title=f'T1<T2={c1}({len(bad)}) T1>T2={c2}({len(good)})')
    fig2.write_html('before_correction_mw_T1.html')
    fig2.show()
    fig3 = px.scatter(x=mw, y=t2, hover_name=t, labels={'x': 'Molecular Weight', 'y': 'T2[s]'}, color=flg,
                      title=f'T1<T2={c1}({len(bad)}) T1>T2={c2}({len(good)})')
    fig3.write_html('before_correction_mw_t2.html')
    fig3.show()


def plot_data2(directory_path):
    data={}
    for filename in os.listdir(directory_path):
        if filename not in ['bmr50553_3.str']: #excluded from output plot
            file_path = os.path.join(directory_path, filename)

            if os.path.isfile(file_path):  # Check if it's a file, not a directory
                data.update(extract_data(file_path))
                print (len(data))
    t1=[]
    dt1=[]
    t2=[]
    dt2=[]
    t=[]
    flg=[]
    good=[]
    bad=[]
    mw=[]
    seq_no=[]
    ent_id=[]
    cc=[]
    unit=[]
    for k in data:
        if 'T1' in data[k] and 'T2' in data[k]:
            if 'T1' in data[k] and 'T2' in data[k]:
                if data[k]['T1'] < data[k]['T2']:
                    c = 'T1<=T2'
                    if k[0] not in bad:
                        bad.append(k[0])
                else:
                    c = 'T2>T2'
                    if k[0] not in good:
                        good.append(k[0])
            unit.append(f'{data[k]['t1unit']}-{data[k]['t2unit']}')
            t1.append(data[k]['T1'])
            mw.append(data[k]['mw'])
            ent_id.append(data[k]['ent_id'])
            seq_no.append(int(k[2]))
            dt1.append(data[k]['T1_err'])
            t2.append(data[k]['T2'])
            dt2.append(data[k]['T2_err'])
            if data[k]['T2']!=0:
                cc.append(data[k]['T1']/data[k]['T2'])
            else:
                cc.append(0.0)
            t.append(k)
            flg.append(c)
    c1=flg.count('T1<=T2')
    c2= len(flg)-c1

    # fig = px.scatter(x=t1,y=t2,hover_name=t,labels = {'x':'T1[s]', 'y' : 'T2[s]'},color=flg,title=f'T1<T2={c1}({len(bad)}) T1>T2={c2}({len(good)})')
    # fig.write_html('after_correction.html')
    # fig.show()
    # fig2 = px.scatter(x=mw,y=t1,hover_name=t,labels = {'x':'Molecular Weight', 'y' : 'T1[s]'},color=flg,title=f'T1<T2={c1}({len(bad)}) T1>T2={c2}({len(good)})')
    # fig2.write_html('after_correction_mw_T1.html')
    # fig2.show()
    # fig3 = px.scatter(x=mw, y=t2, hover_name=t, labels={'x': 'Molecular Weight', 'y': 'T2[s]'}, color=flg,
    #                   title=f'T1<T2={c1}({len(bad)}) T1>T2={c2}({len(good)})')
    # fig3.update_layout(scattermode="group")
    # fig3.write_html('after_correction_mw_t2.html')
    # fig3.show()
    # fig3 = px.scatter(x=seq_no, y=t2,  error_y=dt2, trendline="ewm",trendline_options=dict(halflife=2),hover_name=t, labels={'x': 'Sequence Number', 'y': 'T2[s]'}, color=ent_id,
    #                   title=f'T1<T2={c1}({len(bad)}) T1>T2={c2}({len(good)})')
    # fig3.write_html('after_correction_sq_t2.html')
    # fig3.show()
    # fig3 = px.scatter(x=seq_no, y=t1, error_y=dt1, trendline="ewm", trendline_options=dict(halflife=2),hover_name=t,
    #                   labels={'x': 'Sequence Number', 'y': 'T1[s]'}, color=ent_id,
    #                   title=f'T1<T2={c1}({len(bad)}) T1>T2={c2}({len(good)})')
    # fig3.write_html('after_correction_sq_t1.html')
    # fig3.show()
    # fig3 = px.scatter(x=mw, y=cc, hover_name=t, labels={'x': 'Molecular Weight', 'y': 'T1/T2'}, color=flg,
    #                   title='')
    # fig3.update_layout(scattermode="group")
    # fig3.write_html('after_correction_mw_cc.html')

def extract_data(str_file):
    data = {}
    ent = pynmrstar.Entry.from_file(str_file)
    seq = ent.get_tag('Entity.Polymer_seq_one_letter_code')[0].replace("\n", "")


    print(ent.entry_id, seq)
    mw=calculate_molecular_weight(seq)
    entry_id = ent.entry_id
    for saveframe in ent:
        if saveframe.category in ['heteronucl_T1_relaxation', 'heteronucl_T2_relaxation']:
            t1_unit = saveframe.get_tag('_Heteronucl_T1_list.T1_val_units')
            t2_unit = saveframe.get_tag('_Heteronucl_T2_list.T2_val_units')
            if len(t1_unit) != 0: t1unit = t1_unit[0]
            if len(t2_unit) != 0: t2unit = t2_unit[0]
            print(t1_unit, t2_unit)
            for loop in saveframe:
                if loop.category in ['_T1', '_T2']:
                    for row in loop.data:
                        res = row[loop.tags.index('Comp_ID')]
                        seq = row[loop.tags.index('Comp_index_ID')]
                        atm = row[loop.tags.index('Atom_ID')]

                        if loop.category == '_T1':
                            lst = row[loop.tags.index('Heteronucl_T1_list_ID')]
                            key = (entry_id,lst,seq,res,atm)
                            if key not in data:
                                data[key]={'mw':mw}
                            try:
                                t1=float(row[loop.tags.index('Val')])
                                t1_err = float(row[loop.tags.index('Val_err')])
                                data[key]['T1']=t1
                                data[key]['T1_err']=t1_err
                                data[key]['ent_id']=f'{entry_id}-{lst}'
                                data[key]['t1unit']=t1unit
                            except ValueError:
                                pass
                        else:
                            lst = row[loop.tags.index('Heteronucl_T2_list_ID')]
                            key = (entry_id, lst, seq, res, atm)
                            if key not in data:
                                data[key]={}
                            try:
                                t2=float(row[loop.tags.index('T2_val')])
                                t2_err = float(row[loop.tags.index('T2_val_err')])
                                data[key]['T2'] = t2
                                data[key]['T2_err'] = t2_err
                                data[key]['ent_id'] = f'{entry_id}-{lst}'
                                data[key]['t2unit']=t2unit
                            except ValueError:
                                pass
    return  data



# Protein Molecular Weight Calculator

# Dictionary of average molecular weights of amino acids (in Daltons)
amino_acid_weights = {
    'A': 89.09,  'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
    'E': 147.13, 'Q': 146.15, 'G': 75.07,  'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}

def calculate_molecular_weight(sequence):
    """
    Calculate the molecular weight of a given protein sequence.

    :param sequence: str, protein sequence (single-letter amino acid codes)
    :return: float, molecular weight in Daltons
    """
    weight = sum(amino_acid_weights.get(aa, 0) for aa in sequence)
    return weight





def read_files_in_directory(directory_path):
    """Reads files in a directory one by one."""

    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)

        if os.path.isfile(file_path):  # Check if it's a file, not a directory
           check_entry_for_relaxation_data(file_path)


#read_files_in_directory('./data/input')
#
#plot_data2('./data/output')
plot_data('./data/input')
