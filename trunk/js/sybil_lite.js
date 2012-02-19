// this object is used to maintain the dynamic URL that can be copied by
//  the user as they navigate the interface.
var trackedURL = new Object();
    trackedURL.host = window.location.host;
    trackedURL.pathname = getAppPath();
    
    // these are set later
    trackedURL.project = '';
    trackedURL.orgsOrder = '';
    trackedURL.orgsSelected = '';
    trackedURL.scaffold = '';
    trackedURL.range = '';

function setTrackedURL() {
    trackedURL.project         = $('#project').val();
    trackedURL.orgsOrder       = $('#orgsOrder').val();
    trackedURL.orgsSelected    = $('#orgsSelected').val();
    trackedURL.scaffold        = $('#scaffold').val();
    trackedURL.range           = $('#range').val();
    
    $('#link_view_url').val( 'http://' + trackedURL.host + trackedURL.pathname + 
                             'sybil_lite.cgi?' +
                             'project=' + trackedURL.project + '&' +
                             'orgsOrder=' + trackedURL.orgsOrder + '&' +
                             'orgsSelected=' + trackedURL.orgsSelected + '&' +
                             'scaffold=' + trackedURL.scaffold + '&' +
                             'range=' + trackedURL.range );
    
    //alert("I just set the text.");
}

// gets the path of this application on the server.  I had to use
//  this instead of the window.location properties, which proved
//  to be unreliable.
function getAppPath() {
    var pathArray = location.pathname.split('/');
    var appPath = "/";
    for(var i=1; i<pathArray.length-1; i++) {
        appPath += pathArray[i] + "/";
    }
    return appPath;
}

function checkForm( theForm ) {
    
    storeOrgSpecifications();
    
    var frmStr = $('#frm_panel').serialize();
    
    //alert("form string to submit: " + frmStr);
    
    // clear the result panel and display the loading indicator
    $('#result_panel').empty();
    showProgressIndicator();
    
    $.ajax({
        url: './draw_alignment.cgi',
        dataType: 'json',
        data: frmStr,
        success: function(data, textStatus, jqXHR) {
            //alert("got a successful response back!");
            drawAlignmentsFromJSON(data);
            
            // display the control elements
            $(".alignment_control_elements").show();
            $("#quick_search_c").show();
            $("#btn_save_region_menu").show();
            $("#btn_link_view").show();
            
            // set the range form control value to that returned by the call
            // these retain any scroll/zoom
            $('#range').val( data.range );
            $('#pixelsPerKb').val( data.pixels_per_kb );
            
            setTrackedURL();
        },
        error: function(jqXHR, textStatus, errorThrown){
            alert("Failed to redraw alignment panel! textStatus: (" + textStatus +
                  ") and errorThrown: (" + errorThrown + ")");
                  
            hideProgressIndicator();
        }
    });
}

// uses the current settings of the display panel but shifts to the
//  coordinates of a given gene_id, with some flanking area on the 
//  sides (determined by the current zoom level)
function navigateToGene(gene_id) {
    // get the coordinates of the gene
    //  the AJAX call within this has a callback that does the rest
    $.ajax({
        url: './get_gene_coordinates.cgi',
        dataType: 'json',
        data: { project: $('#project').val(),
                gene: gene_id },
        success: function(data, textStatus, jqXHR) {
            navigateToGene_cb(data);
        },
        error: function(jqXHR, textStatus, errorThrown){
            alert("Failed to get gene coordinates for gene_id: " + gene_id);
        }
    });
}

// callback function that handles most of navigateToGene()
function navigateToGene_cb( gene ) {
    var new_range = { fmin: 0, fmax: 0 };

    // get the coordinates of the existing view (to keep the current scale)
    var range_string = $('#range').val();
        range_string.replace(/ /g,'');
    
    // if range was empty, define some default range at the beginning
    if ( range_string == '' ) {
        range_string = '1-50000';
    }
    
    var old_range_parts = range_string.split(/-/);
    
    // make sure we're dealing with integers here
    old_range_parts[0] = parseInt( old_range_parts[0] );
    old_range_parts[1] = parseInt( old_range_parts[1] );
    gene.fmin = parseInt( gene.fmin );
    gene.fmax = parseInt( gene.fmax );
    
    // calculate the new range
    
    // if the requested gene is larger than the current view, just set our own
    if ( (gene.fmax - gene.fmin) >= old_range_parts[1] - old_range_parts[0] ) {
        new_range.fmin = gene.fmin - 500;
        new_range.fmax = gene.fmax + 500;
        
    } else {
        // else define a range the same size as the current view but with
        //  our gene of interest in the middle of it.
        var flanking_distance = Math.ceil(((old_range_parts[1] - old_range_parts[0]) - 
                                 (gene.fmax - gene.fmin)) / 2);
        
        new_range.fmin = gene.fmin - flanking_distance;
        new_range.fmax = gene.fmax + flanking_distance;
    }
    
    // make sure the range didn't go out of bounds
    if ( new_range.fmin < 0 ) {
        new_range.fmin = 1;
    }
    
    // set the molecule and new range in the form
    $('#range').val( new_range.fmin + '-' + new_range.fmax );
    
    // set the molecule
    $('#scaffold').val( gene.org_abbrev + ';' + gene.molecule );
    
    // submit the form
    checkForm( document.frm_panel );
}

function showProgressIndicator() {
    $('#progress_indicator_panel').show();
}

function hideProgressIndicator() {
    $('#progress_indicator_panel').hide();
}

function drawAlignmentsFromJSON(resp) {

    $('#current_image').val( resp.image_file );
    
    // need to add the image map now - "map_points" json object
    //  each with label => $name, x1 => $x1, y1 => $y1, x2 => $x2, y2 => $y2
    var diagram_html = '<MAP name="featureMap" id="featureMap"></MAP>';
        
    diagram_html += '<img src="show_png.cgi?image=' + resp.image_file + '&keep_image=1" alt="syn_images" usemap="#featureMap">'

    $('#result_panel').append(diagram_html);
    
    var next_area_id_num = 1;
      
    $.each( resp.map_points, function(i, item) {
        var coords_string = item.x1 + ',' + item.y1 + ',' + item.x2 + ',' + item.y2;
    
        $('<area/>', {   "id"       : "area" + next_area_id_num,
                         "shape"    : "rect",
                         "coords"   : coords_string,
                         "title"    : item.label
                    } ).appendTo("#featureMap");
        
        $('#area' + next_area_id_num).hover( 
            function() {
                $('#gene_product').html( item.label );
            }, 
            function() {
                $('#gene_product').html( '' );
            } );
            
        next_area_id_num++;
        
    });
    
    hideProgressIndicator();
}

function storeOrgSpecifications( ) {
    var theForm = document.frm_panel;

    // rebuild the organism order text box (orgsOrder)
    theForm.orgsOrder.value = '';
    
    var organisms = new Array();
    var oselected = new Array();
    
    $('#org_list li').each(
        function() {
            var organism = $(this).text();
                organism = organism.replace(/\s/g,'');
            
            organisms.push( organism );
            
            if ( document.getElementById(organism + '_enabled').value == 1 ) {    
                oselected.push( organism );
            }
        }
    );
    
    theForm.orgsOrder.value = organisms.join(',');
    theForm.orgsSelected.value = oselected.join(',');
}

function addOrg( org_label ) {
    // hide the add button and display the remove button
    document.getElementById(org_label + '_add').style.display = 'none';
    document.getElementById(org_label + '_remove').style.display = 'inline';
    
    document.getElementById(org_label + '_enabled').value = 1;
}

function removeOrg( org_label ) {
    // hide the remove button and display the add button
    document.getElementById(org_label + '_remove').style.display = 'none';
    document.getElementById(org_label + '_add').style.display = 'inline';
    
    document.getElementById(org_label + '_enabled').value = 0;
}

function exportAsSVG() {
    storeOrgSpecifications();
    $('#export_as').val('svg');
    $('#frm_panel').submit();
}

function scrollPanel( direction ) {
    
    $('#scroll_direction').val(direction);

    checkForm( document.frm_panel );
    
    // set the scroll direction back to the default so multiple scrolls don't
    // happen accidentally
    $('#scroll_direction').val('none');
}

function zoomPanel( direction ) {

    $('#zoom_direction').val(direction);

    checkForm( document.frm_panel );
    
    // set the zoom direction back to the default so multiple zooms don't
    // happen accidentally
    $('#zoom_direction').val('none');
    
}

function disableButton( btn ) {
    $(btn).removeClass("ui-state-default");
    $(btn).addClass("ui-state-disabled");
}

function enableButton( btn ) {
    $(btn).removeClass("ui-state-disabled");
    $(btn).addClass("ui-state-default");
}


// saves a region of interest in JSON and PNG files
function saveRoi() {

    storeOrgSpecifications();
    
    var frmStr = $("#frm_panel,#roi_save_label,#roi_save_desc").serialize();

    // load any regions of interest and populate the display
    $.ajax({
        url: './save_roi.cgi',
        dataType: 'json',
        data: frmStr,
        success: function(data, textStatus, jqXHR) {
            saveRoiResponse(data);
        },
        error: function(jqXHR, textStatus, errorThrown){
            alert("Failed to save regions of interest! textStatus: (" + textStatus +
                  ") and errorThrown: (" + errorThrown + ")");
                  
            // hideProgressIndicator();
        }
    });
}

function saveRoiResponse( resp ) {
    $("#save_roi_panel").slideToggle("fast");
    $("#roi_save_label").val('');
    $("#roi_save_desc").val('');
    enableButton( $("#roi_save_cancel") );
    enableButton( $("#roi_save_submit") );
}

// this reads through the regions of interest JSON data returned by the
//  parse_roi.cgi script and transforms it into DOM elements for display
//  within the result panel.  The resultant HTML looks like this:
/*
    <div id="examples_c">
        <div style="background-image: url('./data/Schizos/regions_of_interest/example_inversion_150h.png';);" id="example_inversion" class="alignment_example">
            <h2>Short inversion amid a syntenic region</h2>
            <div class="comment_c">
                <p>
                    This shows a region of good synteny across three genomes but with an 
                    inversion in SO3 starting at coordinate 2315k, made up of 6 genes.  
                    The apparent deletion in SJ1 of the gene found in S03 at 2349K may 
                    just be an assembly artifact where the gene is between existing contigs.
                </p>
            </div>
       </div>
        ... one div per ROI
    </div>
*/

function drawRegionsOfInterestGraphical( resp ) {
    
    $("<div/>", { "id":"examples_c" }).appendTo("#splash_c");

    if ( resp.regions.length > 0 ) {
        $.each(resp.regions, function(key, roi) { 
            $('<div/>', { "class"   : "alignment_example", 
                          "id"      : roi.id
                        } ).appendTo("#examples_c");
            
            $('#' + roi.id).css('background-image', 
                                'url(./data/' + roi.project + '/regions_of_interest/' +
                                 roi.id + '_150h.png)');
            
            $('<h2/>', { text : roi.label } ).appendTo("#" + roi.id);

            $('<div/>', { "class" : "comment_c" } ).appendTo("#" + roi.id);

            $('<p/>', {text : roi.description} ).appendTo("#" + roi.id + " div.comment_c" );

            $("#" + roi.id).click(function() {
                // set the form elements to the data within this ROI
                $('#orgsOrder').val( roi.orgsOrder.join(',') );

                $('#orgsSelected').val( roi.orgsSelected.join(',') );
                for ( var i in roi.orgsSelected ) {
                    addOrg( roi.orgsSelected[i] );
                }

                $('#scaffold').val( roi.scaffold );
                $('#range').val( roi.range );
                $('#flank').val( roi.flank );
                $('#pixelsPerKb').val( roi.pixelsPerKb );

                // submit the alignment request
                checkForm( document.frm_panel );
            });
        });
    } else {
        $('<p/>', {text:"No regions of interest have been defined."}).appendTo("#examples_c");
    }
    
    hideProgressIndicator();
}


// this reads through the regions of interest JSON data returned by the
//  parse_roi.cgi script and transforms it into DOM elements for display
//  within the result panel.  The resultant HTML looks like this:
/*
    <div id="examples_c">
        <table id="examples_c">
            <tr id="roi_row_1">
            
            </tr>
            ... one per ROI
        </table>
    </div>
*/

function drawRegionsOfInterestTabular( resp ) {
    
    $("<div/>", { "id":"examples_c" }).appendTo("#splash_c");
    $("<table/>", { "id":"examples_tbl" }).appendTo("#examples_c");

    if ( resp.regions.length > 0 ) {
        var next_roi_id = 1;
    
        $.each(resp.regions, function(key, roi) { 
            var row_id = 'roi_row_' + next_roi_id++;
            $('<tr/>', { "id" : row_id } ).appendTo("#examples_tbl");
            $('<td/>', { text : roi.label } ).appendTo("#" + row_id);
            $('<td/>', { text : roi.scaffold } ).appendTo("#" + row_id);
            $('<td/>', { text : roi.range } ).appendTo("#" + row_id);
            $('<td/>', { "id" : row_id + '_desc'} ).appendTo("#" + row_id);
            $('<a/>', { "title" : roi.description,
                         text   : '[description]' } ).appendTo( "#" + row_id + '_desc' );
            
            $("#" + row_id).click(function() {
                // set the form elements to the data within this ROI
                $('#orgsOrder').val( roi.orgsOrder.join(',') );

                $('#orgsSelected').val( roi.orgsSelected.join(',') );
                for ( var i in roi.orgsSelected ) {
                    addOrg( roi.orgsSelected[i] );
                }

                $('#scaffold').val( roi.scaffold );
                $('#range').val( roi.range );
                $('#flank').val( roi.flank );
                $('#pixelsPerKb').val( roi.pixelsPerKb );

                // submit the alignment request
                checkForm( document.frm_panel );
            });
        });
    } else {
        $('<p/>', {text:"No regions of interest have been defined."}).appendTo("#examples_c");
    }
    
    hideProgressIndicator();
}







/***************************************************************************
* UseDefaultTextboxValues()
*
* PURPOSE:
* Default Input Values - Handle default values for form elements. Simply call
* function from instead ready scope to automatically handle default values
* for form textbox elements with the ".usehint" class.
*
*   from: http://www.mattvanandel.com/956/jquery-add-smart-textbox-default-values/
**************************************************************************/
function UseDefaultTextboxValues(){
    jQuery('input[type="textbox"].usehint').each(function(){
        //We need to save the element for the sake of the data() method
        var definput = jQuery(this);
        
        //Save the default value for this element
        definput.data('DefVal',definput.val())
        
        //Assign a starting color of silver
        .css('color','silver')
        
        //When the this element gets focus...
        .focus(function(){
            //If the default text is the current value...
            if(definput.val()==definput.data('DefVal')){
                //Clear the text and change color to black
                definput.val('').css('color','black');
            }
        })
    
        //When this element loses focus...
        .blur(function(){
            //If the field is essentially empty...
            if(definput.val().replace(' ','')==''){
                //Use our saved default value
                definput.val(definput.data('DefVal'))

                //And change color back to silver
                .css('color','silver');
            }
        });
        
        //Additional functionality goes here
    });
}



$(function() {
	$("#org_list").sortable();
	$("#org_list").disableSelection();

    $( "#organism_order_d" ).dialog({
        autoOpen: false,
        beforeClose: function(event, ui) { 
                        $('#organism_order_btn').removeClass("ui-state-active"); 
                     }
    });
    $( "#adv_opts_d" ).dialog({
        autoOpen: false,
        beforeClose: function(event, ui) { 
                        $('#adv_opts_btn').removeClass("ui-state-active"); 
                     }
    });
    $( "#molecule_selection_d" ).dialog({
        autoOpen: false,
        beforeClose: function(event, ui) { 
                        $('#molecule_selection_btn').removeClass("ui-state-active"); 
                     }
    });

    UseDefaultTextboxValues();

    // this is a dirty, dirty hack, but JQuery keeps moving my form
    //  dialogs outside of the <form> element, so we're moving
    //  them back here.json format comment
    // http://forum.jquery.com/topic/dialog-takes-fields-out-of-form
    $( "#organism_order_d" ).closest('.ui-dialog').detach().appendTo('form');
    $( "#adv_opts_d" ).closest('.ui-dialog').detach().appendTo('form');
    $( "#molecule_selection_d" ).closest('.ui-dialog').detach().appendTo('form');



    $( "#adv_opts_btn" ).click( function() {
                                        if ( $("#adv_opts_d").dialog('isOpen') ) {
                                            $( "#adv_opts_d" ).dialog('close');
                                        } else {
                                            $( "#adv_opts_d" ).dialog('open');
                                        }
                                     } );

    $( "#molecule_selection_btn" ).click( function() {
                                        if ( $("#molecule_selection_d").dialog('isOpen') ) {
                                            $( "#molecule_selection_d" ).dialog('close');
                                        } else {
                                            $( "#molecule_selection_d" ).dialog('open');
                                        }
                                     } );
    
    $( "#organism_order_btn" ).click( function() {
                                        if ( $("#organism_order_d").dialog('isOpen') ) {
                                            $( "#organism_order_d" ).dialog('close');
                                        } else {
                                            $( "#organism_order_d" ).dialog('open');
                                        }
                                     } );

    $( "#submit_btn" ).click( function() {
                                        checkForm( document.forms["frm_panel"] );
                                     } );

    $( "#btn_save_region_menu" ).click( function () {
                                            $("#save_roi_panel").slideToggle("fast");
                                        });
    $( "#btn_link_view" ).click( function () {
                                            $("#link_view_panel").slideToggle("fast");
                                        });    
    
    $("input#search_term").autocomplete({
        source: "./search_genomes.cgi?project=" + $("#project").val(),
        select: function( event, ui ) {
                    navigateToGene( ui.item.value );
                }
    });
    
    $( "#roi_save_cancel" ).click( function () {
                                $("#roi_save_label").val('');
                                $("#roi_save_desc").val('');
                                $("#save_roi_panel").slideToggle("fast");
                            });
    
    $( "#roi_save_submit" ).click( function () {
                                // don't do anything if it's already disabled
                                if (! $(this).hasClass("ui-state-disabled") ) {
                                    disableButton(this);
                                    disableButton( $("#roi_save_cancel") );
                                    saveRoi();
                                }
                            });
    
    // things with tooltips
    $("#roi_row_1 a[title]").tooltip();
    $("area[title]").tooltip();
    
    // NOTE: later might style control elements with:
    //  http://www.tuttoaster.com/enhancing-forms-using-jquery-ui/3/

    //all hover and click logic for buttons
	$(".fg-button:not(.ui-state-disabled)")
	.hover(
		function(){ 
			$(this).addClass("ui-state-hover"); 
		},
		function(){ 
			$(this).removeClass("ui-state-hover"); 
		}
	)
	.mousedown(function(){
			$(this).parents('.fg-buttonset-single:first').find(".fg-button.ui-state-active").removeClass("ui-state-active");
			if( $(this).is('.ui-state-active.fg-button-toggleable, .fg-buttonset-multi .ui-state-active') ){ $(this).removeClass("ui-state-active"); }
			else { $(this).addClass("ui-state-active"); }	
	})
	.mouseup(function(){
		if(! $(this).is('.fg-button-toggleable, .fg-buttonset-single .fg-button,  .fg-buttonset-multi .fg-button') ){
			$(this).removeClass("ui-state-active");
		}
	});
});


// things within this will start once the DOM is loaded.

// if the URL API is used the form with all its hidden elements is submitted
//  else it displays the regions of interest for the current project
$(document).ready(function() {
    
    if ( $("#url_api_used").val() == 1 ) {
    
        checkForm( document.forms["frm_panel"] );
    
    } else {
    
        // clear the result panel and display the loading indicator
        //$('#result_panel').empty();
        showProgressIndicator();
    
        $("#url_api_used").val('');

        // load any regions of interest and populate the display
        $.ajax({
            url: './parse_roi.cgi',
            dataType: 'json',
            data: { project: $('#project').val() },
            success: function(data, textStatus, jqXHR) {
                if ( $('#graphical_roi_view').val() == 1 ) {
                    drawRegionsOfInterestGraphical(data);
                } else {
                    drawRegionsOfInterestTabular(data);
                }
            },
            error: function(jqXHR, textStatus, errorThrown){
                alert("Failed to draw regions of interest! textStatus: (" + textStatus +
                      ") and errorThrown: (" + errorThrown + ")");

                hideProgressIndicator();
            }
        });
    }
});    




















