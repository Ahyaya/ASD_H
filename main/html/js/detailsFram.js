function toggleTIM(){
    var TIM=document.getElementById("TIM");
    var TIMcheck=document.getElementById("TIMcheck");
    var winList=document.getElementById("secWinForm");
    if(TIMcheck.checked){
        TIM.innerHTML="\u25BC";
        winList.style.display="table";
    }else{
        TIM.innerHTML="\u25BA";
        winList.style.display="none";
    }
}