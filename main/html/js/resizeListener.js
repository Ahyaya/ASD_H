var observer = new MutationObserver(function () {
    var sizeContainer=document.getElementById("sizeContainer");
    if(sizeContainer.offsetHeight > document.documentElement.clientHeight){
        sizeContainer.style.transform="translateX(-50%) translateY(0)";
        sizeContainer.style.top="16px";
    }else{
        sizeContainer.style.transform="translateX(-50%) translateY(-50%)";
        sizeContainer.style.top="50%";
    }
});
var option = {
    childList: true,
    attributes: true,
    subtree: true
}
observer.observe(document.getElementById("sizeTrigger"), option);