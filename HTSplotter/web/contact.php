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

		<title>HTSplotter - Feedback</title>
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
							<a class="nav-link" href="analysis.php">Run analysis</a>
						</li>
						<li class="nav-item">
							<a class="nav-link" href="images/HTSplotterManual.pdf" target="_blank">Manual</a>
						</li>
						<li class="nav-item">
							<a class="nav-link active" href="contact.php">Feedback</a>
						</li>
					</ul>
				</div>
			</nav>
		</header>

		<main role="main" class="container flex-shrink-0">
			<div class="container">
				<div class="row">
					<h2>Contact us</h2>
				</div>
				<div class="row">&nbsp;</div>
				<div class="row">
					<form class="col-12" id="contactform">
						<div class="form-group row">
							<label for="name" class="col-sm-6 col-form-label">Name</label>
							<div class="col-sm-3">
								<input name="fname" id="fname"type="text" class="form-control" placeholder="First name">
							</div>
							<div class="col-sm-3">
								<input name="lname" id="lname"type="text" class="form-control" placeholder="Last name">
							</div>
						</div>
						<div class="form-group row">
							<label for="email" class="col-sm-6 col-form-label">Email address</label>
							<div class="col-sm-6">
								<input type="email" class="form-control" name="email" id="email" aria-describedby="emailHelp">
								<small id="emailHelp" class="form-text text-muted">We'll never share your email with anyone else.</small>
							</div>							
						</div>
						<div class="form-group row">
							<label for="message" class="col-sm-6 col-form-label">Message</label>
							<div class="col-sm-6">
								<textarea class="form-control" name="message" id="message" rows="8"></textarea>
							</div>
						</div>
						<div class="form-group row">
							&nbsp;
						</div>
						<div class="form-group row" id="sendmessage">
							Your message has been sent. Thank you!
						</div>
						<div class="form-group row">
							&nbsp;
						</div>
						<div class="form-group row">
							<div class="col-sm-12 text-right">
                      			<button class="btn btn-primary" id="reset_message">Reset</button>
								<button type="submit" id="send_message" class="btn btn-primary">Submit</button>
							</div>
						</div>
					</form>
				</div>
			</div>
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
