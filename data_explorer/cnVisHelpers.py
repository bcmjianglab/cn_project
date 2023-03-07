def annotatedColor(anObj, df, key):
    '''change color scheme'''
    if key+"_colors" in df.keys():
        print('key already exist')
    else:
        key_index = anObj.obs[key].cat.categories
        key_color = anObj.uns[key+"_colors"]
        kvpair = {k_:v_ for k_, v_ in zip(key_index,key_color)}
        df[key+"_colors"] = df.apply(lambda r: kvpair[r[key]], axis=1)
    return df