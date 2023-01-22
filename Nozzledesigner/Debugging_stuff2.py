def dict_updater(self, eql_nz_props, fzn_nz_props):
    for dictkey, dictval in eql_nz_props.items():
        if dictval != []:
            eqlvalue = dictval[0] * self.prop_conv[dictkey]

            self.eql_nz_props[dictkey].append(eqlvalue)
    for dictkey, dictval in fzn_nz_props.items():
        if dictval != []:
            fznvalue = dictval[0] * self.prop_conv[dictkey]
            self.fzn_nz_props[dictkey].append(fznvalue)