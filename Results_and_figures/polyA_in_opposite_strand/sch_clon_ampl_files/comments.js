//SAmod DIGILOG translit start
var rus_lr2 = ('Е-е-О-о-Ё-Ё-Ё-Ё-Ж-Ж-Ч-Ч-Ш-Ш-Щ-Щ-Ъ-Ь-Э-Э-Ю-Ю-Ю-Ю-Я-Я-Я-Я-Х-Ей-Ей-Ой-Ой-Ий-Ий-Ий-Ий-ё-ё-ж-ч-ш-щ-э-ю-ю-я-я-х-ей-ой-ий-ий').split('-');
var lat_lr2 = ('/E-/e-/O-/o-ЫO-Ыo-ЙO-Йo-ЗH-Зh-ЦH-Цh-СH-Сh-ШH-Шh-ъ'+String.fromCharCode(35)+'-ь'+String.fromCharCode(39)+'-ЙE-Йe-ЙU-Йu-ЫU-Ыu-ЙA-Йa-ЫA-Ыa-Кh-Еy-ЕY-Оy-ОY-Иy-ИY-Ыy-ЫY-ыo-йo-зh-цh-сh-шh-йe-йu-ыu-йa-ыa-кh-еy-оy-иy-ыy').split('-');
var rus_lr1 = ('А-Б-В-Г-Д-Е-З-И-Й-К-Л-М-Н-О-П-Р-С-Т-У-Ф-Х-Х-Ц-Щ-Ы-Я-а-б-в-г-д-е-з-и-й-к-л-м-н-о-п-р-с-т-у-ф-х-х-ц-щ-ъ-ы-ь-ь-я').split('-');
var lat_lr1 = ('A-B-V-G-D-E-Z-I-J-K-L-M-N-O-P-R-S-T-U-F-H-X-C-W-Y-Q-a-b-v-g-d-e-z-i-j-k-l-m-n-o-p-r-s-t-u-f-h-x-c-w-'+String.fromCharCode(35)+'-y-'+String.fromCharCode(39)+'-'+String.fromCharCode(96)+'-q').split('-');
//SAmod DIGILOG translit end


//SAmod DIGILOG translit start
//==========================================
// TRANSLITIRATE (Main)
//==========================================

function translit(all_text)
{
	var obj_ta = document.REPLIER.Post;

	//----------------------------------------
	// It's IE!
	//----------------------------------------
	if ( (all_text=='0') && (ua_vers >= 4) && is_ie && is_win)
	{
        if (obj_ta.isTextEdit)
		{
            obj_ta.focus();
			var sel = document.selection;
			var rng = sel.createRange();
			rng.colapse;
            if((sel.type == "Text" || sel.type == "None") && rng != null)
			{
				rng.text = dotranslate(rng.text);
			}
		}
        else
        {
            obj_ta.value = dotranslate(obj_ta.value);
        }
	}
	//----------------------------------------
	// It's MOZZY!
	//----------------------------------------

	else if ( (all_text=='0') && obj_ta.selectionEnd )
	{
        var ss = obj_ta.selectionStart;
		var st = obj_ta.scrollTop;
		var es = obj_ta.selectionEnd;

		if (es <= 2)
		{
			es = obj_ta.textLength;
		}

		var start  = (obj_ta.value).substring(0, ss);
		var middle = (obj_ta.value).substring(ss, es);
		var end    = (obj_ta.value).substring(es, obj_ta.textLength);

		//-----------------------------------
		// text range?
		//-----------------------------------

		if (obj_ta.selectionEnd - obj_ta.selectionStart > 0)
		{
			middle = dotranslate(middle);
		}

		obj_ta.value = start + middle + end;

		var cpos = ss + (middle.length);

		obj_ta.selectionStart = cpos;
		obj_ta.selectionEnd   = cpos;
		obj_ta.scrollTop      = st;


	}
	//----------------------------------------
	// It's CRAPPY!
	//----------------------------------------
	else
	{
		obj_ta.value = dotranslate(obj_ta.value);
	}

	obj_ta.focus();

	return;
}

//==========================================
// TRANSLITIRATE (String convertion)
//------------------------------------------
// Original code from translit.ru
// by Igor Ilyin (2002-2004)
//==========================================


function dotranslate(text)
{
    text = text.replace(/(^|\s)(www\.\w+[^\s\[\]]+)/ig, "$1[<--$2-->]");
    text = text.replace(/(^|\s)((http|https|news|ftp):\/\/\w+[^\s\[\]]+)/ig, "$1[<--$2-->]");
    
    var txtnew = " ";
    var symb = 0;
    var subsymb = "";
    var trans = 1;
    for (kk=0;kk<text.length;kk++)
    {
        subsymb = text.substr(kk,1);
        if ((subsymb=="[") || (subsymb=="<"))
        {
            trans = 0;
        }
        if ((subsymb=="]") || (subsymb==">"))
        {
            trans = 1;
        }
        if (trans)
        {
            symb = transsymbtocyr(txtnew.substr(txtnew.length-1,1), subsymb);
        }
        else
        {
            symb = txtnew.substr(txtnew.length-1,1) + subsymb;
        }
        txtnew = txtnew.substr(0,txtnew.length-1) + symb;
    }
    txtnew = txtnew.replace(/(\s)Ет/g, "$1Эт");
    txtnew = txtnew.replace(/(\s)ет/g, "$1эт");
    txtnew = txtnew.replace(/(\s):ыес:(\s)/g, "$1:yes:$2");
    txtnew = txtnew.replace(/(\s):но:(\s)/g, "$1:no:$2");
    txtnew = txtnew.replace(/(\s):мол:(\s)/g, "$1:mol:$2");
    txtnew = txtnew.replace(/(\s):лол:(\s)/g, "$1:lol:$2");
    txtnew = txtnew.replace(/(\s):п(\s)/g, "$1:p$2");
    txtnew = txtnew.replace(/(\s):шуффле:(\s)/g, "$1:shuffle:$2");
    txtnew = txtnew.replace(/(\s):о(\s)/g, "$1:o$2");
    txtnew = txtnew.replace(/(\s):еек:(\s)/g, "$1:eek:$2");
    txtnew = txtnew.replace(/(\s):цонфусед:(\s)/g, "$1:confused:$2");
    txtnew = txtnew.replace(/(\s):роллеэс:(\s)/g, "$1:rolleyes:$2");
    txtnew = txtnew.replace(/(\s):щееп:(\s)/g, "$1:weep:$2");
    txtnew = txtnew.replace(/(\s):мад:(\s)/g, "$1:mad:$2");
    txtnew = txtnew.replace(/(\s):цоол:(\s)/g, "$1:cool:$2");
    txtnew = txtnew.replace(/(\s):умник:(\s)/g, "$1:umnik:$2");
    txtnew = txtnew.replace(/(\s):теапот:(\s)/g, "$1:teapot:$2");
    txtnew = txtnew.replace(/(\s):беер:(\s)/g, "$1:beer:$2");
    txtnew = txtnew.replace(/(\s):юмп:(\s)/g, "$1:jump:$2");
    txtnew = txtnew.replace(/(\s):щалл:(\s)/g, "$1:wall:$2");
    txtnew = txtnew.replace(/(\s):мод:(\s)/g, "$1:mod:$2");
    txtnew = txtnew.replace(/(\s):ок:(\s)/g, "$1:ok:$2");
    txtnew = txtnew.replace(/(\s):ноте:(\s)/g, "$1:note:$2");
    txtnew = txtnew.replace(/\[<--/g, "");
    txtnew = txtnew.replace(/-->\]/g, "");

    return txtnew;
}

//==========================================
// TRANSLITIRATE (Symbol convertion)
//------------------------------------------
// Original code from translit.ru
// by Igor Ilyin (2002-2004)
//==========================================

function transsymbtocyr(pretxt,txt)
{
	var doubletxt = pretxt+txt;
	var code = txt.charCodeAt(0);
	if (!(((code>=65) && (code<=123))||(code==35)||(code==39))) return doubletxt;
	var ii;
	for (ii=0; ii<lat_lr2.length; ii++)
	{
		if (lat_lr2[ii]==doubletxt) return rus_lr2[ii];
	}
	for (ii=0; ii<lat_lr1.length; ii++)
	{
		if (lat_lr1[ii]==txt) return pretxt+rus_lr1[ii];
	}
	return doubletxt;
}


//SAmod DIGILOG translit end

function Preview () {
	var earl = "/forums/index.php?act=prev";
	preview = window.open(
		earl,
		"preview",
		"width=550,height=300,toolbar=no,location=no,directories=no,status,menubar=no,scrollbars,resizable,copyhistory=no"
		);
	return true;
} // end Preview




function ShowHide(id1, id2)
{
	if (id1 != '') toggleview(id1);
	if (id2 != '') toggleview(id2);
}
	
function my_getbyid(id)
{
	itm = null;
	if (document.getElementById) itm = document.getElementById(id);
	else if (document.all) itm = document.all[id];
	else if (document.layers) itm = document.layers[id];
	return itm;
}

function toggleview(id)
{
	if ( ! id ) return;
	if ( itm = my_getbyid(id) )
	{
		if (itm.style.display == "none") my_show_div(itm);
		else my_hide_div(itm);
	}
}

function my_hide_div(itm)
{
	if ( ! itm ) return;
	itm.style.display = "none";
}

function my_show_div(itm)
{
	if ( ! itm ) return;
	itm.style.display = "";
}

//==========================================
// pop up window
//==========================================

function PopUp(url, name, width,height,center,resize,scroll,posleft,postop)
{
	showx = "";
	showy = "";
	
	if (posleft != 0) { X = posleft }
	if (postop  != 0) { Y = postop  }
	
	if (!scroll) { scroll = 1 }
	if (!resize) { resize = 1 }
	
	if ((parseInt (navigator.appVersion) >= 4 ) && (center))
	{
		X = (screen.width  - width ) / 2;
		Y = (screen.height - height) / 2;
	}
	
	if ( X > 0 )
	{
		showx = ',left='+X;
	}
	
	if ( Y > 0 )
	{
		showy = ',top='+Y;
	}
	
	if (scroll != 0) { scroll = 1 }
	
	self.name="ori_window";
	var Win = window.open( url, name, 'width='+width+',height='+height+ showx + showy + ',resizable='+resize+',scrollbars='+scroll+',location=no,directories=no,status=no,menubar=no,toolbar=no');
	Win.focus();
	Win.document.close;

}
