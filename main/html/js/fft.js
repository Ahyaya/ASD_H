function formSubmit(){
    var xFile = document.getElementById("xFile");
    var tFile = document.getElementById("tFile");
    var freqFile = document.getElementById("freqFile");
    var winFuncList = document.getElementsByName("winFunc");
    var winFuncCode = 0;
    for(winFuncCode=0;winFuncCode<winFuncList.length;winFuncCode++){
        if(winFuncList[winFuncCode].checked){
            break;
        }
    }
    if(xFile.files.length == 0){
        alert("no singal data selected");
        return;
    }
    if(!document.getElementById("tFileCheckbox").checked && !document.getElementById("udfDt").value){
        document.getElementById("udfDt").focus();
        return;
    }
    if(document.getElementById("tFileCheckbox").checked && tFile.files.length==0){
        alert("no CLKS specified");
        return;
    }
    if(document.getElementById("freqFileCheckbox").checked && freqFile.files.length==0){
        alert("no FREQ specified");
        return;
    }
    document.getElementById("container").innerHTML='';
    var xRequest = new XMLHttpRequest();
    xRequest.onreadystatechange=function(){
        if(this.readyState==4){
            document.getElementById("ResultSpan").style.visibility="visible";
            document.getElementById("sendButton").style.display="inline-block";
            document.getElementById("clearButton").style.display="inline-block";
            document.getElementById("processingSpan").style.display="none";
            console.log("User-define-dt: "+this.getResponseHeader("User-define-dt"));
            console.log("User-define-winfunc: "+this.getResponseHeader("User-define-winfunc"));
            if(this.status==200){
                document.getElementById("downloadLink").innerHTML="Full ASD";
                document.getElementById("filterData").innerHTML="Filtered";
                document.getElementById("downloadLink").href="/processed/"+this.responseText;
                document.getElementById("filterData").href="/processed/"+this.responseText+"_flt.csv";
                plotData("/processed/"+this.responseText+"_flt.csv");
                document.getElementById("downloadLink").style.color="lightgreen";
                document.getElementById("filterData").style.color="cyan";
            }else{
                document.getElementById("downloadLink").innerHTML="Error";
                document.getElementById("filterData").innerHTML="Error";
                document.getElementById("downloadLink").href="#";
                document.getElementById("filterData").href="#";
                document.getElementById("downloadLink").style.color="red";
                document.getElementById("filterData").style.color="red";
            }
        }
    }
    document.getElementById("sendButton").style.display="none";
    document.getElementById("clearButton").style.display="none";
    document.getElementById("processingSpan").style.display="inline-block";

    var userData = new FormData();
    userData.append("xFile",xFile.files[0]);
    if(tFile.files.length>0 && document.getElementById("tFileCheckbox").checked){
        userData.append("tFile",tFile.files[0]);
    }
    if(freqFile.files.length>0 && document.getElementById("freqFileCheckbox").checked){
        userData.append("freqFile",freqFile.files[0]);
    }
    xRequest.open("POST","/fftUpload",true);
    if(!document.getElementById("udfDt").disabled){
        xRequest.setRequestHeader("User-define-dt",document.getElementById("udfDt").value);
    }
    xRequest.setRequestHeader("User-define-winfunc",winFuncList[winFuncCode].value);
    if(!document.getElementById("alpha").disabled){
        xRequest.setRequestHeader("User-define-alpha",document.getElementById("alpha").value);
    }
    xRequest.send(userData);
    return;
}
function fileClear(){
    document.getElementById("xFile").value='';
    document.getElementById("tFile").value='';
    document.getElementById("freqFile").value='';
    document.getElementById("ResultSpan").style.visibility="hidden";
}
function tFileHook(){
    if(document.getElementById("tFileCheckbox").checked){
        document.getElementById("tFile").disabled='';
        document.getElementById("udfDt").disabled='1';
    }else{
        document.getElementById("tFile").disabled='1';
        document.getElementById("udfDt").disabled='';
    }
}
function freqFileHook(){
    if(document.getElementById("freqFileCheckbox").checked){
        document.getElementById("freqFile").disabled='';
    }else{
        document.getElementById("freqFile").disabled='1';
    }
}
async function plotData(fileLocation) {
    await d3.csv(fileLocation,(d)=>{
        return{
            freq: parseFloat(d.freq),
            ASD: parseFloat(d.ASD),
        };
    }).then((data)=>{
        const width = 600;
        const height = 400;
        const padding = { top: 40, right: 96, bottom: 40, left: 96 };
        const gW = width - padding.left;
        const gH = height - padding.top - padding.bottom;

        const svg = d3.select("#container").append("svg").attr("width", width).attr("height", height);
        const graphics = svg.append("g").attr("transform", "translate(48, 48)");
        
        const xScale = d3.scaleLog()
        .domain([d3.min(data, item=>item.freq),d3.max(data, item=>item.freq)]).range([0, gW]);
        const yScale = d3.scaleLog()
        .domain([d3.min(data, item=>item.ASD),d3.max(data, item=>item.ASD)]).range([gH, 0]);

        graphics.append("g").attr("transform", `translate(0, ${gH})`)
        .call(d3.axisBottom(xScale)).attr("stroke", "white");
        graphics.append("g")
        .call(d3.axisLeft(yScale)).attr("stroke", "white");
        graphics.append("g")
        .call(d3.axisTop(xScale)).attr("stroke", "white");
        graphics.append("g").attr("transform", `translate(${gW}, 0)`)
        .call(d3.axisRight(yScale)).attr("stroke", "white");

        const linePath=d3.line().x(d=>{return xScale(d.freq);}).y(d=>{return yScale(d.ASD);});
        graphics.append("path").attr("d",linePath(data)).attr("stroke","cyan").attr("fill","none");
    });
}
function setAlphaByText(){
    var alphaBar=document.getElementById("alpha");
    var alphaText=document.getElementById("alphaInput");
    if(alphaText.value<2){
        alphaText.value=2;
    }
    if(alphaText.value>16){
        alphaText.value=16;
    }
    alphaBar.value=alphaText.value;
}
function setAlphaByRange(){
    document.getElementById("alphaInput").value=document.getElementById("alpha").value;
}
function winReset(){
    document.getElementById("alphaInput").value=10;
    document.getElementById("alpha").value=10;
    document.getElementById("Hanning").checked=1;
    windowChange();
}
function windowChange(){
    var winFuncList = document.getElementsByName("winFunc");
    var winFuncCode = 0;
    for(winFuncCode=0;winFuncCode<winFuncList.length;winFuncCode++){
        if(winFuncList[winFuncCode].checked){
            break;
        }
    }
    switch (winFuncCode) {
        case 5:
        case 6:
        case 9:
        case 10:
            document.getElementById("alphaInput").disabled='';
            document.getElementById("alpha").disabled='';
            break;
        default:
            document.getElementById("alphaInput").disabled="1";
            document.getElementById("alpha").disabled="1";
            break;
    }
}