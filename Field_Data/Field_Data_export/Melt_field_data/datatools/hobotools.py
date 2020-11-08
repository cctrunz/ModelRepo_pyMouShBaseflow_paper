from pandas import read_csv, DataFrame


def read_hobo_csv(csv_file, all_columns=False):
    """
    Reads data from a csv file exported from HOBOware.

    Parameters
    ----------
    csv_file : string
        A string containing the file name of the csv file to be read.
    all_columns : boolean (optional)
        Determines whether to read in all columns or just ones that we
        search for and relabel
        (RH, DewPt, Abs Pres, Temp, Attached, Stopped, Connected, EOF,
        Cond High Range, 
        Cond Low Range). 
        Default = False

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from HOBO csv file.
    """
    skiprows = 1
    index_col = 1
    parse_dates = True
    df = read_csv(csv_file, skiprows=skiprows, index_col=index_col,
                  parse_dates=parse_dates, na_values=['-888.88', '-888.9'])
    # Convert column names into something nicer
    columns = df.columns
    old_columns = columns
    rename_dict = {}
    cond_count = 0
    solar_count = 1
    """
    If you want to search for a term and replace a column name with something different 
    use python tuples instead, example:

    new_names = (('find_this', 'replace_with_this'),) 

    use it like this:

    for old, new in  new_names:
        if old in label:
            new_name = new
            wantcol = True
    """
    new_names = ['RH', 'Gust', 'Wind Speed',
                 'Wind Direction', 'DewPt', 'Abs Pres', 'Rain']
    cols = (old_columns, new_names)

    del df['#']
    for label in columns:
        new_name = label

        if all_columns == False:
            wantcol = False
        else:
            wantcol = True

        for name in new_names:
            if name in label:
                new_name = name
                wantcol = True

        # Account for multiple Solar Radiation Sensors and name them differently.
        # NOTE: current code only allows for up to two solar rad sensors <JZM>
        if 'Solar' in label:
            if solar_count == 1:
                new_name = 'Solar1'
                solar_count = 2
            elif solar_count == 2:
                new_name = 'Solar2'
                solar_count = 3
            else:
                print(">2 Solar Rad Sensors Detecteds")
            wantcol = True

        if 'Temp' in label:
            new_name = 'Temperature'
            wantcol = True

        if 'Detached' in label:
            if all_columns:
                new_name = 'Detached'
        if 'Attached' in label:
            if all_columns:
                new_name = 'Attached'
        if 'Stopped' in label:
            if all_columns:
                new_name = 'Stopped'
        if 'Connected' in label:
            if all_columns:
                new_name = 'Connected'
        if 'End Of File' in label:
            if all_columns:
                new_name = 'EOF'
        if 'Low Range' in label:
            new_name = 'CondLow'
            cond_count += 1
            wantcol = True
        if 'High Range' in label:
            new_name = 'CondHigh'
            cond_count += 1
            wantcol = True

        if wantcol == True:
            rename_dict[label] = new_name
    # If there is only one conductivity column, we'll label it as 'Cond'
    if cond_count == 1:
        old_names = rename_dict.keys()
        for old_name, new_name in rename_dict.iteritems():
            if 'Cond' in new_name:
                cond_key = old_name
        rename_dict[cond_key] = 'Cond'
    df = df.rename(columns=rename_dict)
    # BUG: errors when running, .itervalues <d:2019-04-11 p:2>
#     if all_columns:
    # Trim out unwanted columns
#         s_dict = {}
#         for col in rename_dict.itervalues():
#             s = df[col]
#             s_dict[col] = s
#         df = df[s_dict]

    return df
