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

tutoList.filter(); // Remove all filters
