#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <pthread.h>
#include "../include/ASD_H.h"

#define MAXRECVLEN 4096

const char HTTP200[]="HTTP/1.1 200 OK\r\nContent-Type: text/html; charset=utf-8\r\nConnection: close\r\n";
const char HTTP400[]="HTTP/1.1 400 Bad Request\r\nConnection: close\r\n";

struct conVar_t{
    pthread_mutex_t mutex;
    int clientOn;
} conVar;

void str_xdeplete(char *str) {
    int sp, ep, wrt;
    if(!str[0]){
        return;
    }
    for(ep=strlen(str)-1; ep>0 && (str[ep]=='\n' || str[ep]=='\r' || str[ep]==0x20); --ep);
    for(sp=0; sp<ep && (str[sp]=='\n' || str[sp]=='\r' || str[sp]==0x20); ++sp);
    if(!sp){
        str[ep+1]=0x00;
        return;
    }
    for(wrt=sp;wrt<ep+1;wrt++){
        str[wrt-sp]=str[wrt];
    }
    str[wrt]=0x00;
}

int decode_HTTP_request(char buffer[]){
    int status=0;
    char xFilePath[256]={0}, tFilePath[256]={0}, freqFilePath[256]={0}, dtHeader[128]={0}, winfuncHeader[64]={0}, alphaHeader[64]={0};
    char *sep=buffer, *segm;
	for(segm=strsep(&sep,"\n");segm!=NULL;segm=strsep(&sep,"\n")){
        if(!strncmp(segm,"User-define-alpha:",18)){
            strcpy(alphaHeader,segm+18);
            continue;
        }
        if(!strncmp(segm,"User-define-dt:",15)){
            strcpy(dtHeader,segm+15);
            continue;
        }
        if(!strncmp(segm,"User-define-winfunc:",20)){
            strcpy(winfuncHeader,segm+20);
            continue;
        }
		if(!strncmp(segm,"Content-Disposition: form-data; name=\"",38)){
            segm+=38;
            if(!strncmp(segm,"xFile.path",10)){
                strsep(&sep,"\n");
                if(segm==NULL) break;
                segm=strsep(&sep,"\n");
                str_xdeplete(segm);
                strcpy(xFilePath,segm);
                continue;
            }
            if(!strncmp(segm,"tFile.path",10)){
                strsep(&sep,"\n");
                if(segm==NULL) break;
                segm=strsep(&sep,"\n");
                str_xdeplete(segm);
                strcpy(tFilePath,segm);
                continue;
            }
            if(!strncmp(segm,"freqFile.path",13)){
                strsep(&sep,"\n");
                if(segm==NULL) break;
                segm=strsep(&sep,"\n");
                str_xdeplete(segm);
                strcpy(freqFilePath,segm);
                continue;
            }
        }
    }
    status = (xFilePath[0]!=0) + 2*(tFilePath[0]!=0) + 4*(freqFilePath[0]!=0);
    buffer[0]=0;buffer[256]=0;buffer[512]=0;buffer[768]=0;buffer[1024]=0;buffer[1152]=0;
    strcpy(buffer,xFilePath);
    strcpy(buffer+256,tFilePath);
    strcpy(buffer+512,freqFilePath);
    strcpy(buffer+768,dtHeader);
    strcpy(buffer+1024,winfuncHeader);
    strcpy(buffer+1152,alphaHeader);
    return status;
}

void * client_func(void* arg){
    int fd = *(int*) arg;

	char buffer[4096], msg[4096], outFilePath[256];
	ssize_t numBytesRead;
	int client, status, code, ptr;
    char *xFilePath=buffer, *tFilePath=buffer+256, *freqFilePath=buffer+512, *dtHeader=buffer+768, *winfuncHeader=buffer+1024, *alphaHeader=buffer+1152;
    fftOption_t opts;
    default_fftOption_init(&opts);
    opts.smoothMode=2;
	
	while ((numBytesRead = read(fd, buffer, 4096)) > 0) {
		/*when read from client, handle it with code*/
		buffer[numBytesRead]=0x00;
		status = decode_HTTP_request(buffer);

        sprintf(outFilePath,"%s_ASD.csv",xFilePath);
        code=(status&0x01)?status>>1:-1;
        opts.winOption.winFuncName=winfuncHeader[0]?atoi(winfuncHeader):0;
        opts.dt=dtHeader[0]?atof(dtHeader):opts.dt;
        opts.winOption.alpha=alphaHeader[0]?atof(alphaHeader):opts.winOption.alpha;
        if((opts.dt)<1e-15){
            code=-1;
        }
        fprintf(stderr,"Gate report: code=%d status=%d\n",code,status);
		switch (code)
		{
            case 0:
                status=ASD_H_auto(outFilePath, NULL, xFilePath, NULL, &opts);
                break;
            case 1:
                status=ASD_H_auto(outFilePath, NULL, xFilePath, tFilePath, &opts);
                break;
            case 2:
                status=ASD_H_auto(outFilePath, freqFilePath, xFilePath, NULL, &opts);
                break;
            case 3:
                status=ASD_H_auto(outFilePath, freqFilePath, xFilePath, tFilePath, &opts);
                break;
            case -1:
			default:
				break;
		}

        if(code<0 || status){
            sprintf(msg,"%sContent-Length: 0\r\n\r\n",HTTP400);
            fprintf(stderr,"400 report: code=%d status=%d\n",code,status);
        }else{
            for(ptr=strlen(outFilePath)-1;ptr>0 && outFilePath[ptr]!='/';--ptr);
            ptr=outFilePath[ptr]=='/'?ptr+1:ptr;
            sprintf(msg,"%sUser-define-dt: %e\r\nUser-define-winfunc: %d\r\nUser-define-alpha: %lf\r\nContent-Length: %d\r\n\r\n%s",HTTP200,opts.dt,opts.winOption.winFuncName,opts.winOption.alpha,strlen(outFilePath+ptr),outFilePath+ptr);
        }

		write(fd, msg, strlen(msg));
	}
	
	/*disconnected by client*/
	close(fd);
    return 0x00;
}

void process_connection(int serverfd){
	int recvfd, voidfd[32], fdptr=0;
	pthread_t threadId;
	struct sockaddr_in clientSockAddr;
	socklen_t clientAddrSize;
	clientAddrSize = sizeof(struct sockaddr_in);
	while (1)
	{
		recvfd = accept(serverfd, (struct sockaddr*)&clientSockAddr, &clientAddrSize);
		if (recvfd < 0) {
			continue;
		}

		fdptr &= 0x1F;
		voidfd[fdptr]=recvfd;
		while(pthread_create(&threadId, NULL, client_func, &voidfd[fdptr])){
			sleep(1);
		}
		pthread_detach(threadId);
		fdptr++;
	}
}

void tcpDaemon(int port){
    int listenfd;
    struct sockaddr_in serverSockAddr;
    socklen_t addrlen;
    if ((listenfd = socket(AF_INET, SOCK_STREAM, 0)) < 0){
        fprintf(stderr,"socket() error.\n");
        return;
    }
    bzero(&serverSockAddr, sizeof(serverSockAddr));
    serverSockAddr.sin_family = AF_INET;
    serverSockAddr.sin_port = htons(port);
    //serverSockAddr.sin_addr.s_addr = htonl(INADDR_ANY);
    serverSockAddr.sin_addr.s_addr = inet_addr("127.0.0.1");
    if(bind(listenfd, (struct sockaddr *)&serverSockAddr, sizeof(serverSockAddr)) < 0){
        fprintf(stderr,"bind() error.\n");
        return;
    }
    
    if(listen(listenfd, 16) == -1){
        fprintf(stderr,"listen() error. \n");
        return;
    }

    process_connection(listenfd);

}

int main(int argc, char *argv[])
{
    pthread_mutex_init(&conVar.mutex,NULL);
    tcpDaemon(27016);
    pthread_mutex_destroy(&conVar.mutex);
    return 0;
}
