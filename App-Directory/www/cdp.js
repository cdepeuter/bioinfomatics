$(document).ready(function(){
	$(".dumb").text("leggo");
	clicked = 0;

	$("#go").on("click", function(){
		clicked += 1;
		$("#clicktext").text("clicked:" + clicked);
	})
})