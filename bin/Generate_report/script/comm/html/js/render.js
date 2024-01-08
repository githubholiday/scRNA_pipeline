
class Header extends HTMLElement {
    constructor() {
        super();
        const template = document.createElement('template');
        template.innerHTML = `
            <style>
              .an-header {
                box-sizing: border-box;
                width: 100%;
                height: 70px;
                padding: 0 12px;
                background: linear-gradient(#EEEEEE, #FFFFFF);
                display: flex;
                position: absolute;
                top: 0;
                left: 0;
                right: 0;
                justify-content: space-between;
                border-bottom: 1px solid #eee;
              }
              .an-header .left {
                width: 30%;
                display: flex;
                align-items: center;
            }
            .an-header .left img {
                object-fit: cover;
            }
            .an-header .right {
                width: 30%;
                display: flex;
                flex-direction: column;
                justify-content: center;
                align-items: end;
            }
            </style>
            <div class="an-header">
                <div class="left">
					<img src="./html/img/logo.png" alt="">
                </div>
                <div class="right">
                    <div></div>
                    <div></div>
                </div>
            </div>
        `
        const content = template.content.cloneNode(true);
        this.attachShadow({mode: 'open'}).appendChild(content);
    }
}

class MultiImg extends HTMLElement {
    constructor() {
        super();
        const template = document.createElement('template');
        const templateNoDownload = document.createElement('template');
        template.innerHTML = `
        <style>
        .multiImg {
            margin-top: 20px;
            position: relative;
        }
        
        .multiImg select {
            margin-left: 10px;
            appearance: none;
            width: 220px;
            height: 38px;
            border-radius: 4px;
            padding: 0 24px 0 12px;
            margin-bottom: 20px;
        }
        .drop .tip {
            z-index: 99;
            top: 18px;
            left: 210px;
            position: absolute;
            height: 0;
            width: 0;
            border-top: 4px black solid;
            border-right: 4px solid transparent;
            border-left: 4px solid transparent;
        }
        .multiImg .activeImg {
            display: flex;
            justify-content: center;
            align-items: center;
        }
        .imgwrapper {
            margin: 0 auto;
            overflow: auto;
            width: 600px;
            max-height: 480px;
        }
        .activeImg img {
            width: 100%;
            flex-shrink: 0
        }        
        .title {
            margin-top: 5px;
            text-align: center;
        }
        .downloadWrapper {
            display: flex;
            justify-content: end;
        }
        .download {
            cursor: pointer;
            color: #fff;
            border-radius: 4px;
            display: flex;
            align-items: center;
            display: inline-block;
            text-align: end;
            padding: 0 10px;
            height: 34px;
            line-height: 34px;
            background-color: #191970;
            text-decoration: none;
        }
        </style>
        <div class="multiImg">
            <select name="imgSelect" id="testSelect">
               
            </select>
            <div class="drop">
                <div class="tip"></div>
            </div>
            <div class="imgwrapper">
                <div class="activeImg">
                    <img id="testImg" src="" alt="">
                </div>
            </div>

            <div class="title"></div>

            <div class="downloadWrapper">
                <a class="download">
                    Download
                </a>
            </div>
      </div>
        `

        templateNoDownload.innerHTML = `
        <style>
        .multiImg {
            margin-top: 20px;
            position: relative;
        }
        
        .multiImg select {
            margin-left: 10px;
            appearance: none;
            width: 220px;
            height: 38px;
            border-radius: 4px;
            padding: 0 24px 0 12px;
            margin-bottom: 20px;
        }
        .drop .tip {
            z-index: 99;
            top: 18px;
            left: 210px;
            position: absolute;
            height: 0;
            width: 0;
            border-top: 4px black solid;
            border-right: 4px solid transparent;
            border-left: 4px solid transparent;
        }
        .multiImg .activeImg {
            display: flex;
            justify-content: center;
            align-items: center;
        }
        .imgwrapper {
            margin: 0 auto;
            overflow: auto;
            width: 600px;
            max-height: 480px;
        }
        .activeImg img {
            width: 100%;
            flex-shrink: 0
        }        
        .title {
            margin-top: 5px;
            text-align: center;
        }
        .downloadWrapper {
            display: flex;
            justify-content: end;
        }
        .download {
            cursor: pointer;
            color: #fff;
            border-radius: 4px;
            display: flex;
            align-items: center;
            display: inline-block;
            text-align: end;
            padding: 0 10px;
            height: 34px;
            line-height: 34px;
            background-color: #191970;
            text-decoration: none;
        }
        </style>
        <div class="multiImg">
            <select name="imgSelect" id="testSelect">
               
            </select>
            <div class="drop">
                <div class="tip"></div>
            </div>
            <div class="imgwrapper">
                <div class="activeImg">
                    <img id="testImg" src="" alt="">
                </div>
            </div>

            <div class="title"></div>
      </div>
        `
        const downloadLink = this.getAttribute('download');
        
        const content = downloadLink ? template.content.cloneNode(true) : templateNoDownload.content.cloneNode(true);
        
        const select = content.getElementById('testSelect');
        const title = select.parentElement.querySelector('.title');
        const download = select.parentElement.querySelector('.download') || {};
        const imgs = this.getAttribute('imgs').split(',');
        const titles = this.getAttribute('titles').split(',');
        const index = this.getAttribute('index');
        imgs.map((src, index) => {
            const option = document.createElement('option');
            option.value = src;
            option.innerText = titles[index] ? titles[index] : src;
            select.appendChild(option)    
            return option;
        })
        
        // 初始化
        content.getElementById('testImg').src = imgs[0];
        title.innerText = `${index} ${titles[0]}` || imgs[0];
        download.href = downloadLink ? downloadLink : imgs[0];
        select.addEventListener('change', e => {
            const img = select.parentElement.getElementsByTagName('img')[0];
            if (imgs.indexOf(e.target.value) !== -1) {
                if (titles[imgs.indexOf(e.target.value)]) {
                    title.innerText = `${index} ${titles[imgs.indexOf(e.target.value)]}`
                } else {
                    title.innerText = `${index} ${e.target.value}`
                }
                
            }
            img.src=e.target.value;
            if (!downloadLink) {
                download.href = imgs[0];
            }
        })

        this.attachShadow({mode: 'open'}).appendChild(content);
    }
}

class multiDownload extends HTMLElement {
    constructor() {
        super();
        const template = document.createElement('template');

        template.innerHTML = `
        <style>
        .download {
            cursor: pointer;
            color: #fff;
            border-radius: 4px;
            display: flex;
            align-items: center;
            display: inline-block;
            text-align: end;
            padding: 0 10px;
            height: 34px;
            line-height: 34px;
            background-color: #191970;
            text-decoration: none;
        }
        select {
            margin-left: 10px;
            appearance: none;
            width: 220px;
            height: 38px;
            border-radius: 4px;
            padding: 0 24px 0 12px;
            margin-right: 20px;
        }
        .drop .tip {
            z-index: 99;
            top: 18px;
            right: 125px;
            position: absolute;
            height: 0;
            width: 0;
            border-top: 4px black solid;
            border-right: 4px solid transparent;
            border-left: 4px solid transparent;
        }
        .multiDownload {
            position: relative;
            display: flex;
            justify-content: end;
            align-items: center;
        }
        </style>
        <div class="multiDownload">
            <select id="downloadSelect">
            
            </select>
            <div class="drop">
                <div class="tip"></div>
            </div>
            <a class="download">
                Download
            </a>
        </div>
        `
        const titleAttribute = this.getAttribute('titles');
        const content = template.content.cloneNode(true);
        const downloadLinks = this.getAttribute('download').split(',');
        const select = content.getElementById('downloadSelect');
        const download = select.parentElement.querySelector('.download');
        const titles = typeof titleAttribute === 'string' ? titleAttribute.split(',') : [];
        const regx = /[A-z & .]*/g;

        downloadLinks.map((src, index) => {
            const option = document.createElement('option');
            const match = src.match(regx);
            option.value = src;
            option.innerText = titles[index] ? titles[index] : match ? match[match.length-2] || src : src;
            select.appendChild(option)    
            return option;
        })

        select.addEventListener('change', e => {
            download.href = e.target.value;
        })

        download.href = downloadLinks[0];


        this.attachShadow({mode: 'open'}).appendChild(content);
    }
}

(function () {
    
    window.customElements.define('an-header', Header);
    window.customElements.define('an-multiimg', MultiImg);
    window.customElements.define('an-multidownload', multiDownload);
    const content = document.createElement('div');
    content.classList.add('an-main')
    content.append(...document.body.childNodes)
    document.body.appendChild(content);

    const headerDiv = document.createElement('an-header');
    document.body.appendChild(headerDiv);

    const tocWrapper = document.createElement('div');
    tocWrapper.classList.add('an-wrapper');

    const toc = document.createElement('div');
    toc.id = 'toc'
    tocWrapper.appendChild(toc);
    document.body.appendChild(tocWrapper);
    $('#toc').tocify({
        context: '.an-main',
    });


    $(document).ready( function () {
        $('table').DataTable({
            paging: false,
            filter: false,
            info: false,
			order: [],
        });
    } );
    
})();