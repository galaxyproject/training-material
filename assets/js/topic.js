// setup tutorial list filtering
var options = {
  valueNames: ['tutorial_title']
};
var tutoList = new List('tutorial_list', options)

tutoList.filter(function(item) {
  if (item.values().id > 1) {
     return true;
  } else {
   return false;
  }
}); // Only items with id > 1 are shown in list

tutoList.filter(); // remove filters

// clear search
$('#clear_search').click(function(){
    sb = document.getElementById('tutorial_search');
    sb.value = '';
    sb.dispatchEvent(new KeyboardEvent('keyup'));
});

// search by tutorial tag
$('.tutorial_tag').click(function() {
    sb = document.getElementById('tutorial_search');
    sb.value = this.id;
    sb.dispatchEvent(new KeyboardEvent('keyup'));
});
