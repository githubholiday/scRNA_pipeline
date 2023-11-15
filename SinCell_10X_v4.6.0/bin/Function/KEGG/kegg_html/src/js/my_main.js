/*$(function(){	

	//二维码
	$("#weixin").mouseover(function () {
                $("#wx").css("display","block");
            });
	$("#weixin").mouseleave(function () {
               $("#wx").css("display","none");
            });	
	
	$('#closecode').click(function(){			//关闭二维码
		$('#code').fadeOut('slow');		
		$(this).fadeOut('slow');		
	});	


});
