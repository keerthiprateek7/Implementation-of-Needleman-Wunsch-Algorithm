#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import copy,collections
file=input("Enter the file pdf file you need not require to give .pdb extension to input")
new_file = file + ".pdb"
data = open(new_file,"r")
output = open("sidechain_output.pdb","w")   #now this would open a new file side chain_output so as to display the contents
data.readline()
array=[]
atom_array=[]
amino_count=0
count =1
empty_count=1
dict={}

def for_atoms():
    global amino_count
    global a,b
    for line in data:
        a=''
        b=''
        new_list=line.split()
        if 'ATOM' in line:
            if(new_list[0]=='ATOM'):
                if (new_list[1]=='1'):
                    amino_count=new_list[5]
                if(len(new_list[4])>1):
                    for u in range(0,len(new_list[4])):
                        if (u<1):
                            a=a+(new_list[4][u])
                        else:
                            b=b+(new_list[4][u])
                    new_list[4]=a
                    new_list.insert(5,b)
                    a=''
                    b=''
                if(len(new_list[7])>6):
                    for u in range(0,len(new_list[7])):
                        if (u<6):
                            a=a+(new_list[7][u])
                        else:
                            b=b+(new_list[7][u])
                    new_list[7]=a
                    new_list.insert(8,b)
                    a=''
                    b=''
                if (len(new_list[9]) > 4):
                    for u in range(0,len(new_list[9])):
                        if (u<4):
                            a=a+(new_list[9][u])
                        else:
                            b=b+(new_list[9][u])
                    new_list[9]=a
                    new_list.insert(10,b)
                    a=''
                    b=''
                if (len(new_list[2])>4):
                    for u in range(0,len(new_list[2])):
                        if (u<4):
                            a=a+(new_list[2][u])
                        else:
                            b=b+(new_list[2][u])
                    new_list[2]=a
                    new_list.insert(3,b)
                atom_array.append(new_list)
    atomby_amino(atom_array,amino_count) 
    return

def for_helix():
    data = open(new_file,"r")
    data.readline()
    helix_array=[]
    global counts
    global amino_count
    global a,b
    temp_list={}
    for line in data:
        a=''
        b=''
        new_list=line.split()
        if 'HELIX' in line:
            if (new_list[0]=='HELIX'):
                amino_count=int(new_list[5])
                lastamino_count = int(new_list[8])
                temp_list[amino_count]=lastamino_count
    temp_list=collections.OrderedDict(sorted(temp_list.items()))
    for_helixatoms(temp_list)

def for_helixatoms(temp_list):
    global counts 
    counts=1
    dict_forca={}
    data = open(new_file,"r")
    data.readline()
    helix_array=[]
    for k in temp_list:
        for y in range(k,temp_list[k]+1):
            for line in data:
                a=''
                b=''
                new_list=line.split()
                if 'ATOM' in line:
                    if (new_list[0]=='ATOM'):
                        if(len(new_list[4])>1):
                            for u in range(0,len(new_list[4])):
                                if (u<1):
                                    a=a+(new_list[4][u])
                                else:
                                    b=b+(new_list[4][u])
                            new_list[4]=a
                            new_list.insert(5,b)
                            a=''
                            b=''
                        if(len(new_list[7])>6):
                            for u in range(0,len(new_list[7])):
                                if (u<6):
                                    a=a+(new_list[7][u])
                                else:
                                    b=b+(new_list[7][u])
                            new_list[7]=a
                            new_list.insert(8,b)
                            a=''
                            b=''
                        if (len(new_list[9]) > 4):
                            for u in range(0,len(new_list[9])):
                                if (u<4):
                                    a=a+(new_list[9][u])
                                else:
                                    b=b+(new_list[9][u])
                            new_list[9]=a
                            new_list.insert(10,b)
                            a=''
                            b=''
                        if (len(new_list[2])>4):
                            for u in range(0,len(new_list[2])):
                                if (u<4):
                                    a=a+(new_list[2][u])
                                else:
                                    b=b+(new_list[2][u])
                            new_list[2]=a
                            new_list.insert(3,b)
                        if (int(new_list[5])==y):
                            helix_array.append(new_list)
                        elif (int(new_list[5])>y):
                            break
    for line in helix_array :
        if (line[2]=='CA'):
            dict_forca[counts]=line
            counts=counts+1
    helix_final(dict_forca)

#this function is used to display the helix strand
def helix_final(dict_forca):
    helix_output = open("helix_output.pdb","w")
    x_count=0
    y_count=0
    z_count=0
    rest_total=0
    final_list=copy.deepcopy(dict_forca)
    t=1
    while(t in dict_forca):
        x_count=float(dict_forca[t][6])+float(dict_forca[t+1][6])+float(dict_forca[t+2][6])+float(dict_forca[t+3][6])
        avg_x=round(x_count/4,3)
        y_count=float(dict_forca[t][7])+float(dict_forca[t+1][7])+float(dict_forca[t+2][7])+float(dict_forca[t+3][7])
        avg_y=round(y_count/4,3)
        z_count=float(dict_forca[t][8])+float(dict_forca[t+1][8])+float(dict_forca[t+2][8])+float(dict_forca[t+3][8])
        avg_z=round(z_count/4,3)
        rest_total=float(dict_forca[t][10])+float(dict_forca[t+1][10])+float(dict_forca[t+2][10])+float(dict_forca[t+3][10])
        avg=round(rest_total/4,3)
        final_list[t][1]=t
        final_list[t][2]='S'
        final_list[t][6]=float(avg_x)
        final_list[t][5]=int(dict_forca[t][5])
        final_list[t][7]=float(avg_y)
        final_list[t][8]=float(avg_z)
        final_list[t][9]=float(dict_forca[t][9])
        final_list[t][10]=float(avg)
        t=t+1
        if (len(dict_forca)-t < 3):
            break
    value=len(final_list)-1
    value2=len(final_list)-2
    del final_list[len(final_list)]
    del final_list[value]
    del final_list[value2]
    for r in final_list:
        a= "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(final_list[r][0],final_list[r][1],final_list[r][2],"",final_list[r][3],final_list[r][4],final_list[r][5],"",final_list[r][6],final_list[r][7],final_list[r][8],final_list[r][9],final_list[r][10],final_list[r][11],"")
        data.close()
        helix_output.write(a+'\n')

def atomby_amino(atom_array,amino_count):
    global array
    global count
    global empty_count
    for line in atom_array :
        list2=line
        if (list2[5]==amino_count):
            if (list2[2] not in ('N','C','CA','O','H')):
                array.append(list2)
            else:
                dict[int(amino_count)]=array
                amino_count=str(int(amino_count)+1)
                array=[]
                if (list2[2] not in ('N','C','CA','O','H')):
                    array.append(list2)
    dict[int(amino_count)]=array
    copydict=copy.deepcopy(dict)
    for j in dict:
        length = len(dict[j])
        if (length==0):
            empty_count=empty_count+1
            continue
        x_cord = float(dict[j][0][6])
        y_cord = float(dict[j][0][7])
        z_cord = float(dict[j][0][8])
        last = float(dict[j][0][10])
        for k in range(1,length):
            x_cord = float(dict[j][k][6]) + x_cord 
            y_cord = float(dict[j][k][7]) + y_cord 
            z_cord = float(dict[j][k][8]) + z_cord 
            last = float(dict[j][k][10]) + last 
        averageofx=round(x_cord/length,3)
        averageofy=round(y_cord/length,3)
        averageofz=round(z_cord/length,3)
        averageforlast=round(last/length,3)
        copydict[j][0][1]=count
        copydict[j][0][2]='S'
        copydict[j][0][6]=averageofx
        copydict[j][0][7]=averageofy
        copydict[j][0][8]=averageofz
        copydict[j][0][5]=int(copydict[j][0][5])
        copydict[j][0][9]=float(copydict[j][0][9])
        copydict[j][0][10]=averageforlast
        count=count+1
        r=copydict[j][0]
        a= "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(r[0],r[1],r[2],"",r[3],r[4],r[5],"",r[6],r[7],r[8],r[9],r[10],r[11],"")
        output.write(a+'\n')
        data.close()
    return

def getbeta():
    data = open(new_file,"r")
    data.readline()
    beta_array=[]
    global countz	
    global amino_count
    global a,b
    beta_list={}
    for line in data:
        a=''
        b=''
        new_list=line.split()
        if 'SHEET' in line:
            if (new_list[0]=='SHEET'):
                amino_count=int(new_list[6])
                lastamino_count = int(new_list[9])
                beta_list[amino_count]=lastamino_count
    if(len(beta_list)==0):
        print("This uarticular file does not have any beta strands in it.")
        exit()
    beta_list=collections.OrderedDict(sorted(beta_list.items()))
    for_betaatoms(beta_list)
    

def for_betaatoms(beta_list):
    global countz
    countz=1
    dict_forbeta={}
    data = open(new_file,"r")
    data.readline()
    beta_array=[]
    for line in data:
        a=''
        b=''
        new_list=line.split()
        if 'ATOM' in line:
            if (new_list[0]=='ATOM'):
                if(len(new_list[4])>1):
                    for u in range(0,len(new_list[4])):
                        if (u<1):
                            a=a+(new_list[4][u])
                        else:
                            b=b+(new_list[4][u])
                    new_list[4]=a
                    new_list.insert(5,b)
                    a=''
                    b=''
                if(len(new_list[7])>6):
                    for u in range(0,len(new_list[7])):
                        if (u<6):
                            a=a+(new_list[7][u])
                        else:
                            b=b+(new_list[7][u])
                    new_list[7]=a
                    new_list.insert(8,b)
                    a=''
                    b=''
                if (len(new_list[9]) > 4):
                    for u in range(0,len(new_list[9])):
                        if (u<4):
                            a=a+(new_list[9][u])
                        else:
                            b=b+(new_list[9][u])
                    new_list[9]=a
                    new_list.insert(10,b)
                    a=''
                    b=''
                if (len(new_list[2])>4):
                    for u in range(0,len(new_list[2])):
                        if (u<4):
                            a=a+(new_list[2][u])
                        else:
                            b=b+(new_list[2][u])
                    new_list[2]=a
                    new_list.insert(3,b)
                if (new_list[2] in ('N','CA','C')):
                    beta_array.append(new_list)
    final_betap(beta_array,beta_list)

    
def final_betap(beta_array,beta_list):
    count=1
    finaldict_forbeta={}
    for k in beta_list:
        for e in range(k,beta_list[k]+1):
            for some in beta_array:
                if (int(some[5])==int(e)):
                    finaldict_forbeta[count]=some
                    count=count+1
    finalprocess_forbeta(finaldict_forbeta)

#this function is used to display the beta strands    
def finalprocess_forbeta(finaldict_forbeta):
    beta_output = open("beta_output.pdb","w")
    x_count=0
    y_count=0
    z_count=0
    rest_total=0
    finalbeta_list=copy.deepcopy(finaldict_forbeta)
    t=1
    while(t in finaldict_forbeta):
        x_count=float(finaldict_forbeta[t][6])+float(finaldict_forbeta[t+1][6])+float(finaldict_forbeta[t+2][6])+float(finaldict_forbeta[t+3][6])
        avg_x=round(x_count/4,3)
        y_count=float(finaldict_forbeta[t][7])+float(finaldict_forbeta[t+1][7])+float(finaldict_forbeta[t+2][7])+float(finaldict_forbeta[t+3][7])
        avg_y=round(y_count/4,3)
        z_count=float(finaldict_forbeta[t][8])+float(finaldict_forbeta[t+1][8])+float(finaldict_forbeta[t+2][8])+float(finaldict_forbeta[t+3][8])
        avg_z=round(z_count/4,3)
        rest_total=float(finaldict_forbeta[t][10])+float(finaldict_forbeta[t+1][10])+float(finaldict_forbeta[t+2][10])+float(finaldict_forbeta[t+3][10])
        avg=round(rest_total/4,3)
        finalbeta_list[t][1]=int(t)
        finalbeta_list[t][2]='S'
        finalbeta_list[t][6]=float(avg_x)
        finalbeta_list[t][5]=int(finaldict_forbeta[t][5])
        finalbeta_list[t][7]=float(avg_y)
        finalbeta_list[t][8]=float(avg_z)
        finalbeta_list[t][9]=float(finaldict_forbeta[t][9])
        finalbeta_list[t][10]=float(avg)
        value = t+1
        value1=t+2
        del finalbeta_list[value]
        del finalbeta_list[value1]
        t=t+3
        if (len(finaldict_forbeta)-t < 5):
            break
    value2 = t
    value3=t+1
    value4=t+2
    del finalbeta_list[value2]
    del finalbeta_list[value3]
    del finalbeta_list[value4]
    for r in finalbeta_list:
        a= "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(finalbeta_list[r][0],finalbeta_list[r][1],finalbeta_list[r][2],"",finalbeta_list[r][3],finalbeta_list[r][4],finalbeta_list[r][5],"",finalbeta_list[r][6],finalbeta_list[r][7],finalbeta_list[r][8],finalbeta_list[r][9],finalbeta_list[r][10],finalbeta_list[r][11],"")
        data.close()
        beta_output.write(a+'\n')
for_atoms()
for_helix()
getbeta()

