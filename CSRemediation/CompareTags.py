import pynmrstar

def compar_tags(pdb_nmr_file,bmrb_nmr_file):
    bmrb_tags,bmrb_tags_with_values = get_tags(bmrb_nmr_file)
    pdb_tags, pdb_tags_with_values = get_tags(pdb_nmr_file)
    print (len(bmrb_tags),len(pdb_tags))
    present = []
    missing =[]
    c=0
    x=0
    d=0
    for tag in bmrb_tags:
        if tag in pdb_tags:
            print (f'{",".join(tag)}, {tag in pdb_tags}, {bmrb_tags_with_values[tag] == pdb_tags_with_values[tag]}')
            c+=1
            if bmrb_tags_with_values[tag] != pdb_tags_with_values[tag]: d+=1
        else:
            print(f'{",".join(tag)}, {tag in pdb_tags},None')
            x+=1
    print (len(bmrb_tags),len(pdb_tags),x,c,d)



def get_tags(star_file):
    star = pynmrstar.Entry.from_file(star_file)
    tags=[]
    tags_with_values = {}
    for sf in star:
        for tag in sf.tags:
            tags.append(('saveframe',sf.category,sf.tag_prefix,tag[0]))
            if tags[-1] not in tags_with_values:
                tags_with_values[tags[-1]]=tag[1]
            else:
                #print (f'{tags[-1]} already present in tags list,{sf.name},{tag[1]}',star_file)
                pass
        for lp in sf.loops:
            for t in lp.tags:
                tags.append(('loop',sf.category,lp.category,t))
                if tags[-1] not in tags_with_values:
                    tags_with_values[tags[-1]]=lp.get_tag(t)
    return tags,tags_with_values

if __name__=="__main__":
    compar_tags('../data/30ib_nmr-data.str','../data/bmr35048_3.str')










