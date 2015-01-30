
def quick_dict_to_SVG(Table_dict, width=500, L_width=200, x_location=0, y_location=0):
    title_color = (173,216,230) #lightblue
    header_color = (213,224,197) #a grey green
    font_family="Arial"
    font_size=12 #font_size=font_size
    header_font_size = 14
    line_width=1 #, line_width=line_width
    #good site for color translations: http://www.yellowpipe.com/yis/tools/hex-to-rgb/color-converter.php
    TL = [x_location, y_location]
    row_horz_shift = 5 #shift the row text over to line up with the header text
    display_text = ''
    for header, values in Table_dict.items():
        entries = list(values.keys())
        #Types = values['_type_'] 
        heading = True #default since most row groups have a heading
        #There are at least one entry unter each header
        for entry, explained in values.items():
            if entry =='_type_': 
                #This entry indicates what type of group of rows this is.
                
                #Select the header color
                if explained == 'Heading':
                    bg = header_color
                elif explained == 'Title':
                    bg = title_color
                else: #for comment or other headerless row groups
                    heading = False #
                    
                if heading:    #Build a header for the table
                    head=Table_Header(text=header, width=width, TL=list(TL), background=bg, font_family=font_family, size=header_font_size, line_width=line_width)
                    bg=header_color
                    display_text += head.get_SVG_header() #add the header
                    TL[1]=head.bottom #Move the top left cordinates
                    spacer = Table_rows(font_size=3, width=width, top_left=list(TL), line_width=line_width)
                    spacer.set_count(1)
                    display_text += spacer.get_SVG_rows() #add a spacer after the header
                    TL[1]=spacer.bottom #Move the top left cordinates
                    
            elif "_Comment_starts_at_" in entry: #This is a block of comments
                
                #make a row with one column and multiple lines of text
                txt = clean_text('\n'.join(explained))
                text_list = [[txt]] #
                rows = Table_rows(top_left=list(TL),width=width,font_family=font_family,font_size=font_size, line_width=line_width)
                rows.set_text_list(text_list)
                rows.x_shift=row_horz_shift
                display_text += rows.get_SVG_rows()
                TL[1]=rows.bottom #Move the top left cordinates
                
            elif "_Multiline_Flag_" in entry: #this is an example that spans multiple lines
                
                #Make a row with two columns and multiple lines of text
                multi_left=[]
                multi_right=[]
                #Multiline row entries are lists of tuple pairs
                for line in explained:
                    left = clean_text(line[0])
                    multi_left.append(left)
                    right = clean_text(line[1])
                    multi_right.append(right)
                multi_left_text = '\n'.join(multi_left)
                multi_right_text = '\n'.join(multi_right)
                text_list = [[multi_left_text,multi_right_text]]
                rows = Table_rows(top_left=list(TL),width=width,font_family=font_family,font_size=font_size, line_width=line_width)
                rows.set_text_list(text_list)
                rows.column_locations=[0,L_width]
                rows.x_shift=row_horz_shift
                display_text += rows.get_SVG_rows()
                TL[1]=rows.bottom      #Move the top left cordinates  
                
            else: #this must be rows of examples and explanation
                rows = Table_rows(top_left=list(TL),width=width,font_family=font_family,font_size=font_size, line_width=line_width)
                text_list = [[clean_text(entry), clean_text(explained)]]
                rows.set_text_list(text_list)
                rows.column_locations=[0,L_width]
                rows.x_shift=row_horz_shift
                display_text += rows.get_SVG_rows()
                TL[1]=rows.bottom #Move the top left cordinates
        
        #Add a spacer at the bottom of the row group
        bot_spacer = Table_rows(font_size=3, width=width, top_left=list(TL), line_width=line_width)
        bot_spacer.set_count(1)
        display_text += bot_spacer.get_SVG_rows()
        TL[1]=bot_spacer.bottom #Move the top left cordinates
        
    return (display_text, width+2*line_width+TL[0], TL[1]+2*line_width)