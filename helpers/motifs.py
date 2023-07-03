"""Script handling motifs for evaluation"""

import pandas as pd

class Motif():
    """Motif class contains motif name, its sequence and its id used for indicating it
    """

    def __init__(self, name, seq, id, regex_str=None) -> None:

        # if we pass regex_string, use that instead of str for matching
        if regex_str is None:
            self.regex = seq
        else:
            self.regex = regex_str

        self.seq = seq
        self.name = name
        self.id = id
        self.where = None
        assert len(seq) != 0, "Pass a sequence of length>0"
 
    def __len__(self) -> int:
        # use sequence length
        return len(self.seq)

    def __str__(self):
        return self.name + " " + self.seq


    def ranges(self):
        # all indices
        motif_indices = self.where.nonzero()[0]
        motif_ranges = []
        i = 0

        #print(motif.where)
        #print(motif_indices)
        #print(self)
        #print(len(self))
        #print(motif_indices[-10:])
        while (i < len(motif_indices)):
            # if complete range we take it
            if motif_indices[i + len(self) - 1] - motif_indices[i] == len(self)-1:
                # looks goo, add the motif to list
                motif_ranges.append((int(motif_indices[i]),int(motif_indices[i]+len(self))))
                i += len(self)
            else: 
                print("what")
                print("The given size probably does not match the regex")
                # should not happen if non overlapping motif indications
                # this motif not complete
                i+=1
            
        return motif_ranges

class MotifHandler():
    """Handles all motifs in df
    """

    def __init__(self, motifs) -> None:

        # convert all motif tuples to len 4
        motifs = [m if len(m)==4 else (m[0],m[1],m[2],None) for m in motifs]

        self.df = pd.DataFrame(motifs, columns=["name","seq","id", "regex_str"])
        self.dict = {k:v for (k,v) in zip(list(self.df.seq),list(self.df.id))}
        self.name_seq_dict = {k:v for (k,v) in zip(list(self.df.name),list(self.df.seq))}
        
        self.motifs = []
        # create motif objects
        for idx, row in self.df.iterrows():
            self.motifs.append(Motif(**dict(row)))

        self.df["motif"] = self.motifs

    def __len__(self) -> int:
        return len(self.df)
    
    def get_motif(self, seq : str =None,name : str=None,id : int=None) -> Motif:
        """Get motif by defining either sequence, name or its id. Only pass one of these.

        Args:
            seq (str, optional): motif sequence string. Defaults to None.
            name (str, optional): motif name. Defaults to None.
            id (int, optional): motif id used. Defaults to None.

        Returns:
            Motif: returns defined motif instance
        """
        if seq is not None:
            if seq=="non_motif":
                return "non-motif"
            return self.df[self.df["seq"]==seq].iloc[0]["motif"]
            #Motif(**dict(self.df[self.df["seq"]==seq].iloc[0]))
        if name is not None: 
            if name=="non_motif":
                return "non-motif"
            return self.df[self.df["name"]==name].iloc[0]["motif"]
            #Motif(**dict(self.df[self.df["name"]==name].iloc[0]))
        if id is not None:
            return self.df[self.df["id"]==id].iloc[0]["motif"]
            #Motif(**dict(self.df[self.df["id"]==id].iloc[0]))

    def __iter__(self):
        return iter(self.motifs)
    """
    def __iter__(self):
        self.cidx = 0
        return self

    def __next__(self):
        if self.cidx >= self.__len__():
            raise StopIteration
        motif = Motif(**dict(self.df.iloc[self.cidx]))
        self.cidx += 1
        return motif 
    """

motif_overlap = [
    ("EWSR1","GGGGG"),
    ("FUS", "GGGGG"),
    ("TAF15", "GGGGG"),
    ("HNRNPL", "ACACA"),
    ("PABPN1L", "AAAAA"),
    ("TRA2A", "GAAGA"),
    ("PCBP2", "CCCCC"),
    ("RBFOX2", "GCATG"),
    ("TARDBP", "GTATG"),
    ("HNRNPC", "TTTTT"),
    ("TIA1","TTTTT"),
    ("PTBP3", "TTTCT"),
    ("CELF1", "TATGT"),
    ("FUBP3", "TATAT"),
    ("KHSRP", "TGTAT"),
    ("PUM1", "TGTAT"),
    ("KHDRBS2", "ATAAA")
]

# TODO Anas add all top motifs motifs to here and create motif-handler

# create a motif to id mapping, ids can overlap but proteins should not
motifs_tmp = list(set(map(lambda x: x[1], motif_overlap)))
motifs_tmp = dict(zip(motifs_tmp, range(len(motifs_tmp)))) # (motif, id)

#now add ids to the motif_overlap
motif_overlap = list(map(lambda x: (x[0], x[1], motifs_tmp[x[1]]), motif_overlap))

# MotifHandler takes a list of tuples 
# (protein, motif, id, motif_regex)
motifs = MotifHandler(motif_overlap)