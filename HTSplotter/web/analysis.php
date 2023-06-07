<!doctype html>
<html lang="en">
	<head>
	<!-- Required meta tags -->
		<meta charset="utf-8">
		<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

		<!-- Bootstrap CSS -->
		<link rel="stylesheet" href="bootstrap-4.5.2-dist/css/bootstrap.min.css">
		<!-- own CSS -->
		<link rel="stylesheet" href="css/style.css">

		<title>HTSplotter - run analysis</title>
	</head>
	<body class="d-flex flex-column h-100">
		<header>
			<nav class="navbar navbar-expand-md navbar-dark bg-dark fixed-top">
				<a class="navbar-brand" href="index.php">HTSplotter</a>
				<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarsExampleDefault" aria-controls="navbarsExampleDefault" aria-expanded="false" aria-label="Toggle navigation">
					<span class="navbar-toggler-icon"></span>
				</button>

				<div class="collapse navbar-collapse justify-content-end" id="navbarsExampleDefault">
					<ul class="nav">
						<li class="nav-item">
							<a class="nav-link active" href="analysis.php">Run analysis</a>
						</li>
						<li class="nav-item">
							<a class="nav-link" href="images/HTSplotterManual.pdf" target="_blank">Manual</a>
						</li>
						<li class="nav-item">
							<a class="nav-link" href="contact.php">Feedback</a>
						</li>
					</ul>
				</div>
			</nav>
		</header>

		<main role="main" class="container flex-shrink-0">
			<div class="container">
				<div class="row">
					<h2>Run your own analysis</h2>
				</div>
				<div class="row">&nbsp;</div>
				<div class="modal" tabindex="-1" role="dialog" id="confirm_modal">
					<div class="modal-dialog modal-dialog-centered modal-dialog-scrollable modal-lg" role="document">
						<div class="modal-content">
							<div class="modal-header">
								<h5 class="modal-title">Confirm experiment types</h5>
								<button type="button" class="close" data-dismiss="modal" aria-label="Close">
								<span aria-hidden="true">&times;</span>
								</button>
							</div>
							<div class="modal-body"></div>
							<div class="modal-footer">
								<button type="button" class="btn btn-secondary" data-dismiss="modal" id="cancel_confirm">Close</button>
								<button type="button" class="btn btn-primary" id="confirm_types">Confirm</button>
							</div>
						</div>
					</div>
				</div>
				<div class="row">
					<form class="col-12" id="analysisform">
						<div class="form-group row">
							<label for="csv_files" class="col-sm-6 control-label">Input file(s)</label>
                        	<div class="col-sm-6">
								<span class="btn btn-outline-primary btn-file">
									<input id="csv_files" name="csv_files" type="file" class="file" multiple data-show-upload="true" data-show-caption="false" aria-describedby="inputHelp">
								</span>
								<small id="inputHelp" class="form-text text-muted">Upload one or in case of biological replicates analyses more all TXT-files</small>
                                <small id="inputHelp" class="form-text text-muted">File names without space</small>
							</div>
						</div>
						<div class="form-group row">
							<label for="bio_rep" class="col-sm-6 col-form-label">Biological replicate analysis</label>
							<div class="col-sm-6">
								<select class="form-control" name="bio_rep" id="bio_rep">
									<option>no</option>
									<option>yes</option>
								</select>
							</div>
						</div>
						<div class="form-group row">
							<label for="biorep_fn" class="col-sm-6 col-form-label">Biological replicate desired filename</label>
							<div class="col-sm-6">
								<input type="text" class="form-control" name="biorep_fn" id="biorep_fn">
								<small id="biorepHelp" class="form-text text-muted">Leave empty if no replicate analysis</small>
							</div>
						</div>
						<div class="form-group row">
							<label for="exp_effect" class="col-sm-6 col-form-label">Expected effect</label>
							<div class="col-sm-6">
								<select class="form-control" name="exp_effect" id="exp_effect">
									<option>inhibition</option>
									<option>enhancement</option>
								</select>
							</div>
						</div>
						<div class="form-group row">
							<label for="info_readout" class="col-sm-6 col-form-label">Information readout</label>
							<div class="col-sm-6">
								<input type="text" class="form-control" name="info_readout" id="info_readout">
							</div>
						</div>
						<div class="form-group row">
							<label for="readout_unit" class="col-sm-6 col-form-label">Readout unit</label>
							<div class="col-sm-6">
								<input type="text" class="form-control" name="readout_unit" id="readout_unit">
							</div>
						</div>
						<div class="form-group row">
							<label for="syn_ana" class="col-sm-6 col-form-label">Calculate synergism/antagonism</label>
							<div class="col-sm-6">
								<select class="form-control" name="syn_ana" id="syn_ana">
									<option>Bliss</option>
									<option>HSA</option>
									<option>ZIP</option>
								</select>
							</div>
						</div>
						<div class="form-group row">
							<div class="col-sm-12 hidden">
								<input type="text" class="form-control" name="user_input" id="user_input" value=0>
							</div>
						</div>
						<div class="form-group row" id="loader">
							<div class="col-sm-5">&nbsp;</div>
							<div class="col-sm-2 stage">
								<div class="dot-floating"></div>
							</div>
							<div class="col-sm-5">&nbsp;</div>
						</div>
						<div class="form-group row" id="loader2">
							<div class="col-sm-5">&nbsp;</div>
							<div class="col-sm-2 stage">
								<div class="dot-floating"></div>
							</div>
							<div class="col-sm-5">&nbsp;</div>
						</div>
						<div class="form-group row" id="submitanalysis">
							Your analysis has been submitted. Thank you!
						</div>
						<div class="form-group row" id="analysis_html">
						</div>
						<div class="form-group row">
							&nbsp;
						</div>
						<div class="form-group row">
							<div class="col-sm-12 text-right">
								<button class="btn btn-primary" id="reset_analysis">Reset</button>
								<button type="submit" id="submit_analysis" class="btn btn-primary">Run analysis</button>
							</div>
						</div>
					</form>
				</div>
			</div>
            <div class="row">&nbsp;</div>
            <div class="row">&nbsp;</div>
		</main>

		<footer class="footer mt-auto py-3 flex-shrink-0">	
			 <div class="container">
			 	<div class="row">
			 		<div class="col align-self-center">
			 			<a href="https://www.irc.ugent.be/index.php?id=vermeirssenhome" target="_blank" title="CBIGR"><img src="images/CBIGRlogo2.png" class="rounded footer_img" alt="CBIGR"></a>
						<a href="https://www.spelemanlab.org/" target="_blank" title="SpelemanLab"><img src="images/FSPlogo.png" class="rounded footer_img" alt="SpelemanLab"></a>
						<a href="https://oncornalab.ugent.be" target="_blank" title="OncoRNAlab"><img src="images/oncolab_1.png" class="rounded footer_img" alt="OncoRNAlab"></a>
						<a href="https://www.ugent.be/en" target="_blank" title="UGent"><img src="images/UGent.png" class="rounded footer_img" alt="UGent"></a>
			 		</div>
				</div>
			</div>
		</footer>

		<!-- Optional JavaScript -->
		<!-- jQuery first, then Popper.js, then Bootstrap JS -->
		<script src="https://code.jquery.com/jquery-3.6.0.min.js" integrity="sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=" crossorigin="anonymous"></script>
		<script src="https://cdn.jsdelivr.net/npm/popper.js@1.16.1/dist/umd/popper.min.js" integrity="sha384-9/reFTGAW83EW2RDu2S0VKaIzap3H66lZH81PoYlFhbGU+6BZp6G7niu735Sk7lN" crossorigin="anonymous"></script>
		<script src="bootstrap-4.5.2-dist/js/bootstrap.min.js"></script>
		<script src="js/action.js"></script>
	</body>
</html>
