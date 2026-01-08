from imports import * #Imports

def create_date_code(): #Function to create string that encodes date into saved files
    today=datetime.now() #Pull current date

    if today.month<10: #Processing for single digit months
        month_str='0'+str(today.month)
    else:
        month_str=str(today.month)
    
    if today.day<10: #Processing for single digit days
        day_str='0'+str(today.day)
    else:
        day_str=str(today.day)
    
    today_str=str(today.year)[-2:]+month_str+day_str #Adding year digits and concatenating final string

    return(today_str)