function getCookie(name) {
	var cname = name + "=";
	var dc = document.cookie;

	if (dc.length > 0) {
		begin = dc.indexOf(cname);
		if (begin != -1) {
			begin += cname.length;
			end = dc.indexOf(";", begin);
			if (end == -1)  {
				end = dc.length;
			} // end if
			return (dc.substring(begin, end));
		} // end if
	} // end if
	return null;
}

function help_pop(number)
{
	self.name="ori_window";
	win_op = window.open( '/forums/index.php?act=help&amp;CODE=3&amp;HID=' + number,'HelpWindow','width=500, height=400, resizable=1, scrollbars=1');
	win_op.focus();
	win_op.document.close;
}
