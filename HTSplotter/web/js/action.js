// OWN ACTIONS

$( document ).ready(function() {
	// ANALYSIS FORM ACTIONS
	$("#analysisform").submit(function(event){
		event.preventDefault();
	});

	$('body').on('click', '#submit_analysis', function(){ 
		var email = $('#email').val();
		if (email != '') {
			var data = new FormData();
			var params = $("#analysisform").serializeArray();
			$.each(params, function (key, input) {
				data.append(input.name, input.value);
			});
			var file_data = $('input[name="csv_files"]')[0].files;
			for (var i = 0; i < file_data.length; i++) {
				data.append("csv_file[]", file_data[i]);
			}
			$.ajax({
				type: "POST",
				url: "submit_analysis.php",
				data: data,
				processData: false,
				contentType: false,
				beforeSend: function(){
					$("#loader").css('display', 'flex');
				},
				success: function(msg){
					var obj = jQuery.parseJSON(msg);
					if ('errorfile' in obj){
						$("#submitanalysis").html("Something went wrong. Please click <a href='"+obj.errorfile+"' target='_blank'>here</a> to see the error.");
						$("#submitanalysis").addClass("show");
					} else {
						$('#confirm_modal').modal();
						var exp_types = obj.exp_types.split(",");
						var filenames = obj.filenames.split(",");
						var infopath = obj.infofile.split(",");
						var modal_html = '<div class="container-fluid"><form id="confirmform">';
						$.each(filenames, function( key, value ) {
							modal_html += '<div class="form-group row">';
							modal_html += '<label for="exp_type'+key+'" class="col-md-8 col-form-label">'+filenames[key]+'</label>';
							modal_html += '<div class="col-md-4"><input type="text" class="form-control" id="exp_type'+key+'" name="'+filenames[key]+'" placeholder="Experiment type" value="'+exp_types[key]+'"></div>';
							modal_html += '</div>';
							modal_html += '<div class="form-group row">';
							modal_html += '<div class="col-md-4">Check out the <a href='+infopath[key]+' target="_blank">file details</a></div>';
							modal_html += '</div>';
						});
						modal_html += '</form></div>';
						$('#confirm_modal .modal-body').html(modal_html);
						$('#confirm_modal').modal('show')
						$('body').on('click', '#cancel_confirm', function(){
							$("#submitanalysis").removeClass("show");
							$('#analysisform')[0].reset();
							$("#submit_analysis").removeAttr('disabled','disabled');
						});
						$('body').on('click', '#confirm_types', function(){
							var confirmparams = $("#confirmform").serializeArray();
							var matching_types = true;
							$.each(filenames, function( key, value ) {
								if (confirmparams[key]['name'] == filenames[key]) {
									if (confirmparams[key]['value'] != exp_types[key]){
										matching_types = false;
									}
								} else {
									matching_types = false;
								}
							});
							if (matching_types){
								$('#confirm_modal').modal('hide');
								$("#user_input").val(obj.user_input);
								var data2 = new FormData();
								var params2 = $("#analysisform").serializeArray();
								params2.push({name:'exp_types', value:obj.exp_types});
								$.each(params2, function (key, input) {
									data2.append(input.name, input.value);
								});
								// data2.append('exp_type', exp_type);
								$.ajax({
									type: "POST",
									url: "submit_analysis.php",
									data: data2,
									processData: false,
									contentType: false,
									beforeSend: function(){
										$("#loader2").css('display', 'flex');
									},
									success: function(msg){
										var results = jQuery.parseJSON(msg);
										// var path = results.user_input+"results/";
										// var user_string = results.user_input.split("/")[1];
										var zipfile = results.zip
										$("#submitanalysis").html("Analysis finished successfully. Download your results <a href='"+zipfile+"' target='_blank'>here</a>."); 
										$("#submitanalysis").addClass("show");
									},
									complete:function(data){
										$("#loader2").css('display', 'none');
									}
								});
							} else {
								var infopath = obj.infofile;
								$('#confirm_modal').modal('hide');
								$("#submitanalysis").html("Experiment type(s) do not match. <a href='"+infopath+"' target='_blank'>Click here</a> to see the error.<br>Please fix this issue and resubmit.");
								$("#submitanalysis").addClass("show");
							}
						});
						
					}
					// $("#submitanalysis").addClass("show");
					$("#submit_analysis").attr('disabled','disabled');
				},
				complete:function(data){
					$("#loader").css('display', 'none');
				}
			});
		}
	});

	$('body').on('click', '#reset_analysis', function(){
		$("#submitanalysis").removeClass("show");
		$('#analysisform')[0].reset();
		$("#submit_analysis").removeAttr('disabled','disabled');
	});


	// CONTACT FORM ACTIONS
	$('#contactform').submit(function(event){
		event.preventDefault();
	});

	$('body').on('click', '#send_message', function(){  
		var str = $('#contactform').serialize();
		$.ajax({
			type: "POST",
			url: "contact_response.php",
			data: str,
			success: function(msg){
				$("#sendmessage").addClass("show");
				$("#send_message").attr('disabled','disabled');
				$("#errormessage").ajaxComplete(function(event, request, settings){
					if(msg == 'OK') {
						$("#sendmessage").addClass("show");
					} else {
						$("#sendmessage").removeClass("show");
						result = msg;
					}
					$(this).html(result);
				});
			}
		});
	});

	$('body').on('click', '#reset_message', function(){
		$("#sendmessage").removeClass("show");
		$('#contactform')[0].reset();
		$("#send_message").removeAttr('disabled','disabled');
	});

});
