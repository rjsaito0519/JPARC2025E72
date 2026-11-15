// Optional browser regressions: use an existing Playwright installation.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const fs=require('fs'),assert=require('assert');
const path=require('path');
const html=fs.readFileSync(process.env.WEB_HTML||path.resolve(__dirname,'../calibration/mip_web/mip_web.html'),'utf8').replace('__TOKEN__','fixture');
if(!process.env.PLOTLY_JS)throw Error('Set PLOTLY_JS to an existing plotly.min.js file');
const plotly=fs.readFileSync(process.env.PLOTLY_JS);
const settings={parameters:[10,0,.2],range:[-.5,.5],function:'gaussian',range_mode:'auto'};
const fit={parameters:[10,0,.2],range:[-.5,.5],function:'gaussian',success:true,status:0,covariance:3,ndf:8,curve:[[-.5,0],[0,5],[.5,0]],origin:'fixture'};
async function fixture(browser, staged=false){
 const page=await browser.newPage();
 const state={run:12,particle:'Pi',mode:'t0',fingerprint:'fixture',revision:staged?1:0,staged:staged?{'bh2-3-t0':{settings,result:fit}}:{},job:null,notice:'',gate:{hold_id:'hold',reason:''}};
 let curveResolve,fitResolve,curveReadyResolve,fitReadyResolve;
 const curveReady=new Promise(r=>curveReadyResolve=r),fitReady=new Promise(r=>fitReadyResolve=r);
 const controls={state,polls:0,holdCurve:false,holdFit:false,curveReady,fitReady,releaseCurve:()=>curveResolve?.(),releaseFit:()=>fitResolve?.()};
 await page.route('http://fixture.test/**',async route=>{
  const path=new URL(route.request().url()).pathname;
  if(path==='/'){await route.fulfill({contentType:'text/html',body:html});return;}
  if(path==='/plotly.min.js'){await route.fulfill({contentType:'application/javascript',body:plotly});return;}
  let data;
  if(path==='/state'){controls.polls++;data={...state};}
  else if(path==='/overview')data={...state,pages:1,channels:[3,4].map(seg=>({key:`bh2-${seg}-t0`,label:`BH2 seg${seg} T0`,x:[-1,0,1],counts:[1,5,1],bounds:[-10,10],automatic:fit}))};
  else if(path==='/edit')data={...state,automatic:fit,initial:settings,baseline:settings,automatic_note:'fixture'};
  else if(path==='/curve'){
   if(controls.holdCurve){curveReadyResolve();await new Promise(r=>curveResolve=r);}
   data={curve:[[-.5,99],[0,99],[.5,99]]};
  }else if(path==='/fit'){
   if(controls.holdFit){fitReadyResolve();await new Promise(r=>fitResolve=r);}
   data={...fit,receipt:'receipt',settings};
  }else if(path==='/review')data={...state,rows:[{key:'bh2-3-t0',changes:[['p1',0,.1]]}]};
  else throw Error('Unexpected fixture request '+path);
  try{await route.fulfill({contentType:'application/json',body:JSON.stringify(data)});}catch(e){if(!page.isClosed())throw e;}
 });
 await page.goto('http://fixture.test/?run=12&mode=t0');await page.waitForSelector('.tile button');
 return {page,controls};
}
(async()=>{const browser=await chromium.launch({executablePath:process.env.BROWSER_EXECUTABLE||undefined,headless:true,args:['--no-sandbox']});let failures=[];
 const cases=[
 ['old curve cannot replace a different channel',async()=>{const {page:p,controls:c}=await fixture(browser);try{await p.locator('.tile button').first().click();await p.waitForSelector('#value-1');c.holdCurve=true;await p.locator('#value-1').fill('.1');await c.curveReady;await p.getByRole('button',{name:'← 一覧に戻る',exact:true}).click();await p.getByRole('button',{name:'破棄して戻る',exact:true}).click();await p.locator('.tile button').nth(1).click();await p.waitForFunction(()=>document.getElementById('channel-label').textContent==='BH2 seg4 T0');c.releaseCurve();await p.waitForTimeout(200);assert(!(await p.evaluate(()=>document.getElementById('plot').data[1].y)).includes(99),'old channel curve appeared on new channel');}finally{c.releaseCurve();await p.close();}}],
 ['fit prevents editing inputs while result is pending',async()=>{const {page:p,controls:c}=await fixture(browser);try{await p.locator('.tile button').first().click();await p.waitForSelector('#value-1');c.holdFit=true;await p.getByRole('button',{name:'Fit',exact:true}).click();await c.fitReady;assert(await p.locator('#value-1').isDisabled(),'parameter stayed editable during fit');assert(await p.locator('#function').isDisabled(),'method stayed editable during fit');c.releaseFit();await p.waitForFunction(()=>document.getElementById('fit-status').textContent.includes('Fit成功'));assert(await p.locator('#value-1').isEnabled(),'parameter never re-enabled');}finally{c.releaseFit();await p.close();}}],
 ['stale confirmation disables Write after another tab changes staging',async()=>{const {page:p,controls:c}=await fixture(browser,true);try{await p.getByRole('button',{name:'Writeする内容を確認',exact:true}).click();assert(await p.locator('#write').isEnabled());c.state.revision=2;await p.waitForTimeout(2800);assert(await p.locator('#write').isDisabled(),'stale review retained enabled Write');}finally{await p.close();}}]
 ];
 for(const [name,test] of cases){try{await test();console.log('PASS',name);}catch(e){failures.push(name+': '+e.message);console.log('FAIL',name,e.message);}}
 await browser.close();if(failures.length)process.exitCode=1;
})().catch(e=>{console.error(e);process.exit(1)});
