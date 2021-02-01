#Fixes a latex .bib file to display only the specified number of authors.
# Returns an updated latex file that should be usable by the latex file.

from cantrips import recursiveReplaceMultipleSpaces

def makeNewAuthorList(raw_author_entry, max_n_authors, n_authors_before_etAl, abrev_first, abrev_middle):
    raw_author_entry = raw_author_entry.strip()
    print 'raw_author_entry = ' + str(raw_author_entry)  
    raw_author_set = raw_author_entry.split(' and ') #parse the raw author string by the string 'and'
    print 'raw_author_set = ' + str(raw_author_set) 
    new_author_set = []
    for author in raw_author_set:
        author = author.strip() 
        print 'author = ' + author
        if ',' in author:
            #Input author formatted like 'Last, first second'
            last_name = author.split(',')[0].strip() #from start to comma
            first_name = author.split(',')[1].strip().split(' ')[0].strip() #from comma to space
            if len(author.split(',')[1].strip().split(' ')) > 1:
                mid_name = author.split(',')[1].strip().split(' ')[1].strip()  #from space to end
            elif '.' in first_name[0:-1]:
                mid_name = first_name.split('.')[1] + '.'
                first_name = first_name.split('.')[1] + '.'
            else:
                mid_name = ''
        else:
            #input author formatted like 'first second last'
            if len(author.strip().split(' ')) > 2: 
                last_name = author.split(' ')[2].strip() 
                first_name = author.split(' ')[0].strip()
                mid_name = author.split(' ')[1].strip()
            else:
                last_name = author.split(' ')[1].strip()
                first_name = author.split(' ')[0].strip()
                
                mid_name = ''
        print 'first_name = ' + str(first_name)
        print 'last_name = ' + str(last_name)
        print 'mid_name = ' + str(mid_name) 
        if abrev_first:
            first_name = first_name[0].strip() + '.'
        if abrev_middle and len(mid_name) > 0:
            mid_name = mid_name[0].strip() + '.'
        new_author = last_name + ', ' + first_name + ' ' + mid_name
        new_author_set = new_author_set + [new_author]

    print 'new_author_set = ' + str(new_author_set) 
    if len(new_author_set) <= max_n_authors:
        print 'here1'
        authors_to_join = new_author_set
    else:
        print 'here2' 
        authors_to_join = new_author_set[0:n_authors_before_etAl]#  + [ 'et. al.']

    print 'authors_to_join = ' + str(authors_to_join) 
    new_author_entry = ' and '.join(authors_to_join)
    if len(new_author_set) <= max_n_authors: new_author_entry = new_author_entry + ' and et al.'
    
    return new_author_entry 
    
    
        

#Expects every entry to be formatted as a dictionary 
def trimAuthorList(rawLatexFile, updatedLatexFile, max_n_authors = 3, n_authors_before_etAl = 2, abrev_first = 1, abrev_middle = 1):
    in_file = open (rawLatexFile, 'r')
    data = in_file.readlines()
    lines = [line.rstrip('\n') for line in data]
    updated_entries = []
    line_in_entry = 0
    new_entry = []
    n_braces = 0
    in_author = 0
    author_entry = ''
    new_author_entries = []
    new_lines = [] 
    for line in lines:
        if '@' in line:
            #Start of a new actual entry
            print 'starting entry' 
            line_in_entry = 1
        if line_in_entry:
            n_braces = n_braces + len([char for char in line if char == '{']) - len([char for char in line if char == '}'])
            # check if line is something that needs to be updated
            if 'author' in line:
                print line
                in_author = 1
            if in_author:
                author_entry = author_entry + line + ' ' 
                if line.strip()[-2:] == '},' or line.strip()[-1] == '}':
                    #End of author entry
                    author_entry = recursiveReplaceMultipleSpaces(author_entry) 
                    #print 'author_entry0 = ' + str(author_entry)
                    in_author = 0
                    author_entry = author_entry.strip()
                    #print 'author_entry1 = ' + str(author_entry)
                    author_entry = ''.join(author_entry.split('{')[1:])
                    #print 'author_entry2 = ' + str(author_entry)
                    author_entry = ''.join(author_entry.split('}')[0:-1])
                    #print 'author_entry3 = ' + str(author_entry)
                    new_author_entry = '  author = {' + makeNewAuthorList(author_entry, max_n_authors, n_authors_before_etAl, abrev_first, abrev_middle) + '},'
                    print 'new_author_entry = ' + str(new_author_entry)
                    new_author_entries = new_author_entries + [[new_author_entry]]
                    new_lines = new_lines + [new_author_entry]
            else:
                new_lines = new_lines + [line]
            if n_braces == 0:
                print 'Finishing entry'
                line_in_entry = 0
                new_entry = []

    new_file = open(updatedLatexFile, 'w')
    for new_line in new_lines:
        new_file.write("%s\n" % new_line)

    return new_lines 
