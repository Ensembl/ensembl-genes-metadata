process SET_DATE {
    output:
    stdout

    script:
    if(params.date && !params.full_screen)
        """
            echo ${params.date}
        """
    else if(params.full_screen && !params.date)
        """
            ${params.metadata_db.host} ${params.metadata_db.database}  -NB -e "SELECT DATE_FORMAT(date_value, '%m/%d/%Y') from update_date WHERE update_type = 'full_screen';"
        """
    else if(!params.date && !params.full_screen)
        """
            ${params.metadata_db.host} ${params.metadata_db.database}  -NB -e "SELECT DATE_FORMAT(DATE_SUB(date_value, INTERVAL 1 DAY), '%m/%d/%Y') from update_date WHERE update_type = 'regular_update';"
        """ 
    else
        error "Invalid parameters to set up date"

}